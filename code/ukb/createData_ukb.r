suppressPackageStartupMessages(source("wis37138/msm-timeScales-ieb-ic-ckd/code/ukb/helpers_ukb.r"))

# set up paths ----
dir_pheno <- "/wis37138/_transfer/UKBB_level3/"
dir_geno <- "/wis37138/_transfer/ukb594/"
dir_geno2 <- "/wis37138/_transfer/UKBB_eMR_AC/"
dir_rrt <- "/wis37138/_transfer/UKBB_eMR_AC/severe_kidney_events/"
dir_out <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb/"

file_pheno <- "20250605_UKBB_level3_data.txt"
file_geno <- "eGFR_meta_ea.full_cc.addOA.594_allchr_switched2crossSectionalLoweringAllele.doset"
file_geno_12 <- "01_UKBB_eMR_AC_genetics_12decline_SNPs_final_switched2crossSectionalLoweringAllele.doset"
file_geno_rs2075570 <- "02_UKBB_genetics_rs2075570.txt"
file_PCs <- "all_pops_non_eur_pruned_within_pop_pc_covs_id20272_id23940.tsv"
file_coefs <- "/wis37138/kidney_function_decline/results/UKBB/new_n_348275_m_1520382/UKBB_595_summary_7approaches.xlsx" # for PGS

# set flags ----
min_trajectory_length <- 2
exclude_aki_6months <- TRUE
exclude_nephrectomy_6months <- TRUE
exclude_nonGPclinical <- TRUE
death_info <- TRUE
# cutoffs <- c(90, 60, 30)
cutoffs <- c(60, 30, 15)
mid <- TRUE # default: TRUE
imputation <- TRUE # default: TRUE
rounding <- 0 # default: NULL
age_scale <- TRUE # default: TRUE
censoring <- 25 # default: NULL; only for time scale

file_pheno_beforeFlags <- paste0("pheno_ukb_beforeFlags_", paste(cutoffs, collapse = ""), ".txt")
file_pheno_out <- paste0("pheno_ukb_", paste(cutoffs, collapse = ""), ".txt")
file_events_out <- paste0("events_ukb_", paste(cutoffs, collapse = ""), ".txt")

# read data ----

## phenotype data ----
pheno_id <- "id20272"
pheno <- read_file(file.path(dir_pheno, file_pheno)) %>%
    rename(id := !!pheno_id) %>%
    select(-c(record_number, time_diff)) %>%
    distinct()

## gp clinical idx ----
gp_idx <- readRDS(file.path(dir_out, "gp_idx.rds"))

## geno ----
geno_id <- "id23940"
geno <- read_file(file.path(dir_geno, file_geno)) %>%
    rename(!!geno_id := id) %>%
    left_join(
        read_file(file.path(dir_geno2, file_geno_12)) %>% select(one_of(geno_id, "rs28857283_A_G")),
        by = geno_id) %>%
    left_join(
        read_file(file.path(dir_geno2, file_geno_rs2075570)) %>% select(-ID),
        by = geno_id)

colnames(geno)[colnames(geno) == "9:140085182_A_G_A_G"] <- "rs72763290_A_G"

## event data ----
aki <- extract_date(pheno, "AKI")
nephrectomy <- extract_date(pheno, "nephrectomy")
dialysis <- extract_date(pheno, "dialysis")
transplant <- extract_date(pheno, "transplantation")
eskd_code <- extract_date(pheno, "ESKD")
pregnancy <- extract_date(pheno, "pregnancy")
death <- extract_date(pheno, "emr_death")

## risk factors ----
diabetes <- extract_date(pheno, "Diabetes") # not time-varying (single onset date), therefore keep incorporating like this
risk_factors <- c("sc_BMI", "emr_BMI", "sc_UACR", "emr_UACR")
diagnostic_codes <- c("Albuminuria_NOS", "Microalbuminuria", "non_smoking", "smoking")
# pheno_tvc <- create_tvc(pheno, risk_factors, diagnostic_codes)
# pheno_tvc <- pheno_tvc %>%
#   mutate(
#     smoking = case_when(
#       smoking == 1           ~ 1,
#       non_smoking == 1       ~ 0,
#       TRUE                   ~ NA_real_
#     )
#   ) %>%
#   select(-non_smoking)
# fwrite(pheno_tvc, file.path(dir_out, "pheno_tvc_ukb.txt"), sep = "\t", quote = FALSE)
pheno_tvc <- fread(file.path(dir_out, "pheno_tvc_ukb.txt"))

pheno <- pheno %>%
    filter(record_type %in% c("emr_crea", "sc_crea") & !is.na(sex) & !is.na(eGFRcrea) & !is.na(age_at_exam) & !is.na(date_at_exam)) %>%
    mutate(
        sex01 = ifelse(sex == "Female", 0, 1)
    ) %>%
    select(-one_of("record_value", "diagnostic_code", "date_at_measurement", "year_at_measurement", "year_at_exam", "record_type", "sc_bl_or_fu", "time_since_first", "death", "sex01")) %>%
    relocate(eGFRcrea, .after = sex)

eGFR_last <- create_eGFR_final_cutoff(pheno, cutoff = cutoffs[3])

# preprocessing ----

## merge geno and PCs ----
pheno <- pheno %>%
    # include only non-related EUR
    filter((!panrelated) & (panpop == "EUR")) %>% #n=351785; m=1553727
    # merge geno data
    left_join(geno, by = geno_id) %>%
    # exclude individuals without genetic data
    filter(!!as.symbol(geno_id) %in% unlist(geno[geno_id])) %>%
    # exclude measurements taken at age<18
    filter(age_at_exam >= 18)

## merge event data ----

### merge data with a single event per individual
pheno <- pheno %>%
    left_join(
        diabetes %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_diabetes = date_at_exam),
        by = c("id")) %>%
    left_join(
        aki %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_aki = date_at_exam),
        by = c("id")) %>%
    left_join(
        dialysis %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_dialysis = date_at_exam),
        by = c("id")) %>%
    left_join(
        nephrectomy %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_nephrectomy = date_at_exam),
        by = c("id")) %>%
    left_join(
        transplant %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_transplant = date_at_exam),
        by = c("id")) %>%
    left_join(
        eskd_code %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_eskd_code = date_at_exam),
        by = c("id")) %>%
    left_join(
        eGFR_last %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_last_eGFR = date_at_exam),
        by = "id") %>%
    left_join(
        death %>%
            group_by(id) %>% arrange(date_at_exam) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(date_at_exam_death = date_at_exam),
        by = c("id"))

## create diabetes variable (here: no info -> no diabetes) ----
setDT(pheno)
pheno[, date_at_exam := as.IDate(date_at_exam, format="%Y-%m-%d")]
pheno[, diabetes := ifelse(is.na(date_at_exam_diabetes) | date_at_exam < date_at_exam_diabetes, 0, 1)]

## create ESKD and cut trajectories afterwards ----
pheno[, date_at_exam_eskd := pmin(date_at_exam_dialysis, date_at_exam_transplant, date_at_exam_eskd_code, date_at_exam_last_eGFR, na.rm = TRUE)]
pheno <- pheno[, .SD[date_at_exam <= date_at_exam_eskd | is.na(date_at_exam_eskd)], by = id]

## cut trajectories after death ----
pheno <- pheno[, .SD[date_at_exam <= date_at_exam_death | is.na(date_at_exam_death)], by = id]
fwrite(pheno, file.path(dir_out, file_pheno_beforeFlags), sep = "\t", quote = FALSE)

# apply flags ----
# pheno <- fread(file.path(dir_out, file_pheno_beforeFlags))

## exclude assessments +- 6 months around AKI
if(exclude_aki_6months) {
    pheno <- pheno[is.na(date_at_exam_aki) | !((as.integer(date_at_exam_aki) - 180) <= as.integer(date_at_exam) &
                                                as.integer(date_at_exam) <= (as.integer(date_at_exam_aki) + 180))]
    file_pheno_out <- gsub(".txt", "_noAKI6months.txt", file_pheno_out)
    file_events_out <- gsub(".txt", "_noAKI6months.txt", file_events_out)
}

## exclude assessments +- 6 months around nephrectomy
if(exclude_nephrectomy_6months) {
    pheno <- pheno[is.na(date_at_exam_nephrectomy) | !((as.integer(date_at_exam_nephrectomy) - 180) <= as.integer(date_at_exam) &
                                                as.integer(date_at_exam) <= (as.integer(date_at_exam_nephrectomy) + 180))]
    file_pheno_out <- gsub(".txt", "_noNephrectomy6months.txt", file_pheno_out)
    file_events_out <- gsub(".txt", "_noNephrectomy6months.txt", file_events_out)
}

## keep only individuals with at least min_trajectory_length eGFR assessments ----
pheno <- pheno[, if (.N >= min_trajectory_length) .SD, by = id]
file_pheno_out <- gsub(".txt", paste0("_ni", min_trajectory_length,  "+.txt"), file_pheno_out)
file_events_out <- gsub(".txt", paste0("_ni", min_trajectory_length,  "+.txt"), file_events_out)

## exclude individuals who are not in GP clinical ----
if(exclude_nonGPclinical) {
    pheno <- pheno[id %in% gp_idx]
    file_pheno_out <- gsub(".txt", "_GPclinical.txt", file_pheno_out)
    file_events_out <- gsub(".txt", "_GPclinical.txt", file_events_out)
}

# create status, time, etc. ----
pheno <- add_status(pheno, cutoffs)
pheno <- add_vars(pheno)
fwrite(pheno, file.path(dir_out, file_pheno_out), sep = "\t", quote = FALSE)

# create event datasets ----
file_pheno_out <- "pheno_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
pheno <- fread(file.path(dir_out, file_pheno_out))

adjust_ids <- pheno[need_adjustment == TRUE, unique(id)]
age_eskd_dt <- pheno[!is.na(age_eskd), .(age_eskd = unique(age_eskd)), by = id]
# rs_vector <- c(
#   "rs77924615", "rs13334589", "rs434215", "rs28857283", "rs4930319",
#   "rs10224002", "rs1458038", "rs74209810", "rs2783971", "rs854922",
#   "rs2076668", "rs2921093", "rs1047891", "rs35969577", "rs79760705",
#   "rs3812036", "rs2279463", "rs13230509", "rs10851885", "rs9905761",
#   "rs11657044", "rs34446110", "rs35709439", "rs2075570"
# )
# rs_vector <- read.xlsx(file_coefs) %>%
#     filter(Type %in% c("decline", "stable")) %>%
#     pull(SNP)
# pattern <- paste0("^(", paste(rs_vector, collapse = "|"), ")")
# snp_cols <- grep(pattern, colnames(pheno), value = TRUE)
snp_cols <- grep("^rs[0-9]", colnames(pheno), value = TRUE)
selected_cols <- c(
  "id", "agestart", "agestop", "from", "to", "transition", "date_at_exam",
  "status", "episode", "sex", "eGFRcrea", "diabetes", "age0", "agestop_raw",
  snp_cols
)

events_raw <- create_events(pheno, snp_cols = snp_cols, selected_cols = selected_cols, adjust_ids = adjust_ids, age_eskd_dt = age_eskd_dt, death_info = death_info)
saveRDS(events_raw, file.path(dir_out, gsub(".txt", "_beforePGS.rds", file_events_out)))
file_events_out <- "events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
events_raw <- readRDS(file.path(dir_out, file_events_out))

# add PGS ----
weights_decline <- read.xlsx(file_coefs) %>%
    filter(Type =="decline") %>%
    mutate(
      rsid = paste0(SNP, "_", OA, "_", EA),
      weights = interaction_coef_age_rs/sum(interaction_coef_age_rs)*12) %>%
    select(rsid, weights)

weights_umod <- read.xlsx(file_coefs) %>%
    filter(SNP %in% c("rs77924615", "rs13334589", "rs74209810")) %>%
    mutate(
      rsid = paste0(SNP, "_", OA, "_", EA),
      weights = interaction_coef_age_rs/sum(interaction_coef_age_rs)*3) %>%
    select(rsid, weights)

weights_stable <- read.xlsx(file_coefs) %>%
    filter(Type =="stable") %>%
    mutate(
      rsid = paste0(SNP, "_", OA, "_", EA),
      weights = main_coef_age_rs/sum(main_coef_age_rs)*11) %>%
    select(rsid, weights)

weights_cross_595 <- read.xlsx(file_coefs) %>%
    mutate(
      rsid = paste0(SNP, "_", OA, "_", EA),
      weights = cross_coef/sum(cross_coef)*595) %>%
    select(rsid, weights)

weights_cross_594_umod <- read.xlsx(file_coefs) %>%
    filter(SNP != "rs77924615") %>%
    mutate(
      rsid = paste0(SNP, "_", OA, "_", EA),
      weights = cross_coef/sum(cross_coef)*594) %>%
    select(rsid, weights)

# events_raw <- events_raw %>%
#   rowwise() %>%
#   mutate(
#     pgs_decline = sum(c_across(all_of(weights_decline$rsid)) * weights_decline$weights),
#     pgs_umod = sum(c_across(all_of(weights_umod$rsid)) * weights_umod$weights),
#     pgs_stable = sum(c_across(all_of(weights_stable$rsid)) * weights_stable$weights),
#     pgs_cross_595 = sum(c_across(all_of(weights_cross_595$rsid)) * weights_cross_595$weights),
#     pgs_cross_594_umod = sum(c_across(all_of(weights_cross_594_umod$rsid)) * weights_cross_594_umod$weights)
#     ) %>%
#   ungroup()

weight_list <- list(
  decline        = weights_decline,
  umod           = weights_umod,
  stable         = weights_stable,
  cross_595      = weights_cross_595,
  cross_594_umod = weights_cross_594_umod
)

for (nm in names(weight_list)) {
  wt <- weight_list[[nm]]
  events_raw[
    , paste0("pgs_", nm) := rowSums(.SD * wt$weights)
    , .SDcols = wt$rsid
  ]
}

rs_vector <- read.xlsx(file_coefs) %>%
    filter(Type %in% c("decline", "stable")) %>%
    pull(SNP)
pattern <- paste0("^(", paste(rs_vector, collapse = "|"), ")")
snp_cols_selected <- grep(pattern, colnames(pheno), value = TRUE)
non_snp_cols <- setdiff(colnames(events_raw), snp_cols)
cols_to_keep <- c(non_snp_cols, snp_cols_selected)
events_raw <- events_raw[, ..cols_to_keep]

saveRDS(events_raw, file.path(dir_out, gsub(".txt", ".rds", file_events_out)))

## adjustment for interval-censoring ----
# file_events_out <- "events_ukb_906030_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
file_events_out <- "events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
events_raw <- readRDS(file.path(dir_out, gsub(".txt", ".rds", file_events_out)))

if(mid) {
    events <- adjust_events(events_raw, pheno)
    file_events_out <- gsub(".txt", "_mid.txt", file_events_out)
} else {
    events <- events_raw
    file_events_out <- gsub(".txt", "_end.txt", file_events_out)
}

## imputation (for jumps 0->2, 0->3, and 1->2)
if(imputation) {
    events <- impute_events(events, pheno, cutoffs)
    if(sum(events$transition %in% c("0->2", "0->3", "1->3")) > 0) {
      events <- impute_events(events, pheno, cutoffs)
    }
    file_events_out <- gsub(".txt", "_imputed.txt", file_events_out)
}

## rounding (to reduce the number of distinct event times)
if(!is.null(rounding)) {
    if(!age_scale) {
        events <- round_events(events, r = rounding, age_scale = FALSE)
        file_events_out <- gsub(".txt", paste0("_r", rounding, "_timeScale.txt"), file_events_out)
    } else {
        events <- round_events(events, r = rounding, age_scale = TRUE)
        file_events_out <- gsub(".txt", paste0("_r", rounding, "_ageScale.txt"), file_events_out)
    }
}

# ## administrative censoring
# if(!is.null(censoring) & !age_scale) {
#     events <- censor_events(events, c=censoring)
#     file_events_out <- gsub(".txt", paste0("_cens", censoring, ".txt"), file_events_out)
# }

# ## clock-reset
# if(age_scale) {
#     events <- reset_clock(events)
#     file_events_out <- gsub(".txt", "_clockReset.txt", file_events_out)
# }

## saving
fwrite(events, file.path(dir_out, file_events_out), sep = "\t", quote = FALSE)
saveRDS(events, file.path(dir_out, gsub(".txt", ".rds", file_events_out)))
# file_events_out <- "events_ukb_906030_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.txt"
# events <- readRDS(file.path(dir_out, gsub(".txt", ".rds", file_events_out)))

# create ped ----
events_ct <- events %>%
    rename(tstart = agestart, tstop = agestop) %>%
    add_counterfactual_transitions()

ped_events <- as_ped(
  data         = events_ct,
  formula      = Surv(tstart, tstop, status) ~ .,
#   cut          = seq(19, 80, 1),
  transition   = "transition",
  id           = "id",
  censor_code  = 0,
  timescale    = "calendar")

ped_events <- ped_events %>%
  merge_tvc(
    pheno_tvc %>% mutate(age = round_up_half(age, 0)),
    vars = c("sc_BMI", "emr_BMI", "sc_UACR", "Albuminuria_NOS", "Microalbuminuria", "smoking"))

ped_events <- add_timeScales(ped_events)
file_ped_events_out <- paste0("ped_", file_events_out)
saveRDS(ped_events, file.path(dir_out, gsub(".txt", ".rds", file_ped_events_out)))

# # split ped by transition ----
# ped_events_list <- split_ped(ped_events)

# for(i in 1:3) {
#     saveRDS(ped_events_list[[i]], file.path(dir_out, gsub(".txt", paste0("_", i-1, i, ".rds"), file_ped_events_out)))
# }

# create anonymized data ----
set.seed(123)

n <- 20 # number of IDs to sample per strata

# 1) compute counts per id
id_counts <- events %>%
  count(id, name = "n_rows")

# 2) sample 20 ids from each strata
ids_1 <- id_counts %>%
  filter(n_rows == 1) %>%
  sample_n(n) %>%
  pull(id)

ids_2 <- id_counts %>%
  filter(n_rows == 2) %>%
  sample_n(n) %>%
  pull(id)

ids_3 <- id_counts %>%
  filter(n_rows == 3) %>%
  sample_n(n) %>%
  pull(id)

# 3) combine into one vector
sample_n_ids <- c(ids_1, ids_2, ids_3)

events_ct_anon <- events %>%
  # 0) select only the IDs you want to anonymize
  filter(id %in% sample_n_ids) %>%
  # 1) grab just the columns you want
  select(id, agestart, agestop, from, to, transition, status, sex) %>%
  # 2) rename for PED consistency, stash the original ID away as orig_id
  rename(
    tstart  = agestart,
    tstop   = agestop,
    orig_id = id
  ) %>%
  # 3) do your per-ID noise injection
  group_by(orig_id) %>%
  arrange(tstart, .by_group = TRUE) %>%
  group_modify(~ {
    df <- .
    n   <- nrow(df)
    # discretized Gaussian shift on first tstart
    noise_start <- round(rnorm(1, mean = 0, sd = 5))
    # discrete lengths 1–20
    noise_stops <- sample(1:20, n, replace = TRUE)

    new_tstart <- numeric(n)
    new_tstop  <- numeric(n)

    # first interval
    new_tstart[1] <- df$tstart[1] + noise_start
    new_tstop[1]  <- new_tstart[1] + noise_stops[1]

    # chain subsequent intervals
    if (n >= 2) {
      new_tstart[2] <- new_tstop[1]
      new_tstop[2]  <- new_tstart[2] + noise_stops[2]
    }
    if (n >= 3) {
      new_tstart[3] <- new_tstop[2]
      new_tstop[3]  <- new_tstart[3] + noise_stops[3]
    }

    df %>% mutate(
      tstart = new_tstart,
      tstop  = new_tstop
    )
  }) %>%
  ungroup() %>%
  # 4) remap every orig_id → 1,2,3… and re‐randomise sex
  mutate(
    id  = as.integer(factor(orig_id, levels = unique(orig_id))),
    sex = ifelse(rbinom(n(), 1, 0.5) == 1, "Male", "Female")
  ) %>%
  # 5) explicitly select id back into the final tibble
  select(
    id,       # your new anon ID
    tstart,
    tstop,
    from,
    to,
    transition,
    status,
    sex
  ) %>%
  # 6) finally add your counterfactual transitions
  add_counterfactual_transitions()

ped_events_anon <- as_ped(
  data         = events_ct_anon,
  formula      = Surv(tstart, tstop, status) ~ .,
  cut          = seq(min(events_ct_anon$tstart), max(events_ct_anon$tstop), 1),
  transition   = "transition",
  id           = "id",
  censor_code  = 0,
  timescale    = "calendar")

ped_events_anon <- add_timeScales(ped_events_anon)
file_ped_events_anon_out <- "ped_events_anon.rds"
saveRDS(ped_events_anon, file.path(dir_out, file_ped_events_anon_out))
file_events_anon_out <- "events_anon.rds"
saveRDS(events_ct_anon, file.path(dir_out, file_events_anon_out))

a <- readRDS(file.path(dir_out, file_ped_events_anon_out))
b <- readRDS(file.path(dir_out, file_events_anon_out))


# compute interval length per state ----
library(dplyr)
library(lubridate)


pheno <- pheno %>%
  mutate(date_at_exam = as.Date(date_at_exam))


avg_gap_years <- function(df) {
  df %>%
    arrange(id, date_at_exam) %>%
    group_by(id) %>%
    mutate(gap_days = as.numeric(difftime(lead(date_at_exam), date_at_exam, units = "days"))) %>%
    summarise(id_avg_years = mean(gap_days / 365.25, na.rm = TRUE), .groups = "drop") %>%
    filter(!is.nan(id_avg_years)) %>%
    summarise(mean_years = mean(id_avg_years),
              n_ids = n())
}

res_00 <- pheno %>%
  filter(flag1 == 0, flag2 == 0) %>%
  avg_gap_years() %>%
  mutate(scenario = "flag1=0 & flag2=0")

res_10 <- pheno %>%
  filter(flag1 == 1, flag2 == 0) %>%
  avg_gap_years() %>%
  mutate(scenario = "flag1=1 & flag2=0")

res_01plus <- pheno %>%
  filter(flag2 == 1) %>%
  avg_gap_years() %>%
  mutate(scenario = "flag2=1")

summary_gaps <- bind_rows(res_00, res_10, res_01plus) %>%
  select(scenario, mean_years, n_ids)

summary_gaps
