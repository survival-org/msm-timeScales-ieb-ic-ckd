
source("wis37138/kidney_survival/code/helpers.r")

# set up paths ----
dir_pheno <- "/wis37138/_transfer/UKBB_eMR_AC/"
dir_geno <- "/wis37138/_transfer/ukb594/"
dir_linker <- "/wis37138/_transfer/UKBB_eMR_AC/LinkerFile/"
dir_PCs <- "/wis37138/_transfer/UKBB_eMR_AC/"
dir_rrt <- "/wis37138/_transfer/UKBB_eMR_AC/severe_kidney_events/"
dir_out <- "/wis37138/kidney_survival/data/ukb/"

file_pheno <- "01_UKBB_eMR_AC_20230915.txt"
file_pheno_beforeFlags <- "pheno_ukb_beforeFlags.txt"
file_geno <- "eGFR_meta_ea.full_cc.addOA.594_allchr.doset"
file_linker <- "link_id23940_id20272_id8343.txt"
file_PCs <- "all_pops_non_eur_pruned_within_pop_pc_covs_id20272_id23940.tsv"
file_aki <- "03_AKI_id20272.txt"
file_transplant <- "05_kidneyTransplant_id20272.txt"
file_dialysis <- "07_dialysis_id20272.txt"
file_nephrectomy <- "08_Nephrectomy_id20272.txt"
file_eskd <- "06_ESKD_id20272.txt"
file_pregnancy <- "09_Pregnancy_id20272.txt"
file_pheno_out <- "pheno_ukb.txt"
file_events_out <- "events_ukb.txt"
file_events_withBackTransitions_out <- "events_withBackTransitions_ukb.txt"

# set flags ----
min_trajectory_length <- 10
exclude_aki_6months <- TRUE

# read data ----

## phenotype data ----
pheno_id <- "id20272"
pheno <- read_file(file.path(dir_pheno, file_pheno)) %>% # n=454907; m=2102174
    rename(id := !!pheno_id)

## geno ----
geno_id <- "id23940"
geno <- read_file(file.path(dir_geno, file_geno)) %>%
    rename(!!geno_id := id)

## linker for geno id ----
linker <- read_file(file.path(dir_linker, file_linker)) %>%
    select(one_of(pheno_id, geno_id))

## PCs ----
PCs <- read_file(file.path(dir_PCs, file_PCs)) %>%
    distinct(get(pheno_id), .keep_all = TRUE) %>%
    select(one_of(pheno_id), starts_with("PC"))

## event data ----
aki <- read_file(file.path(dir_rrt, file_aki)) %>%
    mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
    filter(!is.na(id20272) & !is.na(event_dt)) %>%
    select(id20272, event_dt)

kidneyTransplant <- read_file(file.path(dir_rrt, file_transplant)) %>%
    mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
    filter(!is.na(id20272) & !is.na(event_dt)) %>%
    select(id20272, event_dt)

dialysis <- read_file(file.path(dir_rrt, file_dialysis)) %>%
    mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
    filter(!is.na(id20272) & !is.na(event_dt)) %>%
    select(id20272, event_dt)

nephrectomy <- read_file(file.path(dir_rrt, file_nephrectomy)) %>%
    mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
    filter(!is.na(id20272) & !is.na(event_dt)) %>%
    select(id20272, event_dt)

eskd_code <- read_file(file.path(dir_rrt, file_eskd)) %>%
    mutate(event_dt = as.Date(event_dt, format = "%d/%m/%Y")) %>%
    filter(!is.na(id20272) & !is.na(event_dt)) %>%
    select(id20272, event_dt)

eGFR15 <- pheno %>%
    group_by(id) %>%
    filter(any(eGFRcrea < 15)) %>%  # Keep only individuals with any eGFRcrea < 15
    filter(eGFRcrea < 15) %>%        # Filter rows with eGFRcrea < 15
    slice_min(order_by = event_dt, n = 1, with_ties = FALSE) %>%  # Get the first date when eGFRcrea < 15
    select(id, event_dt) %>%    # Keep only id and event_dt columns
    ungroup()

# preprocessing ----

## merge geno and PCs ----
pheno <- pheno %>%
    # # include only non-related EUR
    # filter((!panrelated) & (panpop == "EUR")) %>% #n=351785; m=1553727
    # merge geno id
    left_join(linker, by = c(id = pheno_id)) %>%
    # merge PCs
    left_join(PCs, by = c(id = pheno_id))
    # # exclude individuals without genetic data
    # filter(!!as.symbol(geno_id) %in% unlist(geno[geno_id]))

## merge event data ----
pheno <- pheno %>%
    left_join(
        aki %>%
            group_by(!!as.symbol(pheno_id)) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_aki = event_dt),
        by = c(id = "id20272")) %>%
    left_join(
        dialysis %>%
            group_by(!!as.symbol(pheno_id)) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_dialysis = event_dt),
        by = c(id = "id20272")) %>%
    left_join(
        nephrectomy %>%
            group_by(!!as.symbol(pheno_id)) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_nephrectomy = event_dt),
        by = c(id = "id20272")) %>%
    left_join(
        kidneyTransplant %>%
            group_by(!!as.symbol(pheno_id)) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_transplant = event_dt),
        by = c(id = "id20272")) %>%
    left_join(
        eskd_code %>%
            group_by(!!as.symbol(pheno_id)) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_eskd_code = event_dt),
        by = c(id = "id20272")) %>%
    left_join(
        eGFR15 %>%
            group_by(id) %>% arrange(event_dt) %>% slice_head(n = 1) %>% ungroup() %>%
            rename(event_dt_15 = event_dt),
        by = "id")

## create ESKD and cut trajectories afterwards ----
setDT(pheno)
pheno[, event_dt := as.IDate(event_dt, format="%Y-%m-%d")]
pheno[, event_dt_eskd := pmin(event_dt_dialysis, event_dt_transplant, event_dt_eskd_code, event_dt_15, na.rm = TRUE)]
pheno <- pheno[, .SD[event_dt <= event_dt_eskd | is.na(event_dt_eskd)], by = id]
# fwrite(pheno, file.path(dir_out, file_pheno_beforeFlags), sep = "\t", quote = FALSE)

# apply flags ----
pheno <- fread(file.path(dir_out, file_pheno_beforeFlags))

## exclude assessments +- 6 months around AKI
if(exclude_aki_6months) {
    pheno <- pheno[is.na(event_dt_aki) | !(event_dt_aki - 180 <= event_dt & event_dt <= event_dt_aki + 180)]
    file_pheno_out <- gsub(".txt", "_noAKI6months.txt", file_pheno_out)
    file_events_out <- gsub(".txt", "_noAKI6months.txt", file_events_out)
    file_events_withBackTransitions_out <- gsub(".txt", "_noAKI6months.txt", file_events_withBackTransitions_out)
}

## keep only individuals with at least min_trajectory_length eGFR assessments ----
pheno <- pheno[, if (.N >= min_trajectory_length) .SD, by = id]
file_pheno_out <- gsub(".txt", paste0("_ni", min_trajectory_length,  "+.txt"),file_pheno_out)
file_events_out <- gsub(".txt", paste0("_ni", min_trajectory_length,  "+.txt"),file_events_out)
file_events_withBackTransitions_out <- gsub(".txt", paste0("_ni", min_trajectory_length,  "+.txt"),file_events_withBackTransitions_out)

# create time & age0 ----
pheno <- pheno[order(id, event_dt)]
pheno[, t := (as.numeric(event_dt) - min(as.numeric(event_dt))) / as.numeric(365.25), by = id]
pheno[, age0 := age[which.min(t)], by = id]

# create status variable ----
pheno[, flag60 := as.integer(eGFRcrea < 60)]
pheno[, flag30 := as.integer(eGFRcrea < 30)]
pheno[, status := 0]
pheno[event_dt == event_dt_eskd, status := 3]

pheno[, status_withBackTransitions := {
  s <- status  # Initialize status vector for the group

  # Handle eGFRcrea < 30 for at least 2 consecutive measurements (status = 2)
  for (i in seq_len(.N)) {
    if (flag30[i] == 1 && i > 1 && flag30[i - 1] == 1) {
      # Backfill status to the start of the consecutive sequence
      start_i <- i - 1
      while (start_i >= 1 && flag30[start_i] == 1) {
        if (s[start_i] < 2) s[start_i] <- 2
        start_i <- start_i - 1
      }
      if (s[i] < 2) s[i] <- 2
    }
  }

  # Handle eGFRcrea < 60 for at least 2 consecutive measurements (status = 1)
  for (i in seq_len(.N)) {
    if (flag60[i] == 1 && s[i] == 0 && i > 1 && flag60[i - 1] == 1) {
      # Backfill status to the start of the consecutive sequence
      start_i <- i - 1
      while (start_i >= 1 && flag60[start_i] == 1) {
        if (s[start_i] < 1) s[start_i] <- 1
        start_i <- start_i - 1
      }
      if (s[i] < 1) s[i] <- 1
    }
  }

  s  # Return the updated status vector
}, by = id]

pheno[, status := {
    s <- status  # Initialize status vector with current statuses (0 or 3)

    # Assign status 2 where eGFRcrea < 30 for at least 2 consecutive measurements
    idx2 <- which(flag30 == 1 & shift(flag30, type = "lag") == 1)
    s[idx2] <- pmax(s[idx2], 2)
    s[idx2 - 1] <- pmax(s[idx2 - 1], 2)

    # Assign status 1 where eGFRcrea < 60 for at least 2 consecutive measurements
    idx1 <- which(flag60 == 1 & shift(flag60, type = "lag") == 1)
    s[idx1] <- pmax(s[idx1], 1)
    s[idx1 - 1] <- pmax(s[idx1 - 1], 1)

    # Ensure that status does not decrease over time (excluding transitions to status 3)
    for (i in 2:.N) {
        if (s[i - 1] == 3) {
            s[i] <- 3  # Once in status 3, remain in status 3
        } else if (s[i] != 3) {
            s[i] <- max(s[i - 1], s[i])  # Enforce non-decreasing status
        }
    }
    s  # Return the updated status vector
}, by = id]

pheno[, t_eskd := as.numeric(event_dt_eskd - min(event_dt))/as.numeric(365.25), by = id]
pheno[, run_id := rleid(status), by = id]
pheno[, run_id_withBackTransitions := rleid(status_withBackTransitions), by = id]
pheno[, t_last := max(t), by = id]
pheno[, event_dt_last := max(event_dt), by = id]
pheno[, need_adjustment := !is.na(event_dt_eskd) & (event_dt_eskd != event_dt_last), by = id]

fwrite(pheno, file.path(dir_out, file_pheno_out), sep = "\t", quote = FALSE)
# pheno <- fread(file.path(dir_out, file_pheno_out))

# create event datasets ----
adjust_ids <- pheno[need_adjustment == TRUE, unique(id)]
t_eskd_dt <- pheno[!is.na(t_eskd), .(t_eskd = unique(t_eskd)), by = id]

## events without back transitions ----
events <- pheno[, .(
  from = as.character(first(status)),
  tstart = first(t),
  tstop = last(t),
  age0 = first(age0), # time-constant
  sex = first(sex), # time-constant
  # TBD: genetics (all 13+11)
), by = .(id, run_id)]
events[, to := as.character(shift(from, type = "lead")), by = id]
events[, is_last_event := seq_len(.N) == .N, by = id]
events[, need_adjustment := id %in% adjust_ids]
events[is_last_event == TRUE, to := ifelse(
  need_adjustment,
  3,
  ifelse(from %in% c(0, 1, 2), "cens", NA)
)]
events[, transition := paste0(from, "->", to)]
events <- events[!(from == 3 & is.na(to))]
events <- merge(events, t_eskd_dt, by = 'id', all.x = TRUE)
events[need_adjustment == TRUE & is_last_event == TRUE, stop := t_eskd]
events[, episode := seq_len(.N), by = .(id, transition)]
events[, status := ifelse(to == "cens", 0, 1)]
events[, tstop2 := ifelse(
    to == "3",
    t_eskd,
    shift(tstart, type = "lead")
), by = id]
events[is.na(tstop2), tstop2 := tstop]
# TBD
## time since entering state (speak to Johannes + Andreas)
## diabetes (time-varying)
## eGFR at entry into state != 3 (time-varying)

setorder(events, id, tstart)
events <- events[, .(id, tstart, tstop, tstop2, from, to, transition, status, episode)]

fwrite(events, file.path(dir_out, file_events_out), sep = "\t", quote = FALSE)

## events with back transitions ----
events_withBackTransitions <- pheno[, .(
  from = as.character(first(status_withBackTransitions)),
  tstart = first(t),
  tstop = last(t)
), by = .(id, run_id_withBackTransitions)]
events_withBackTransitions[, to := as.character(shift(from, type = "lead")), by = id]
events_withBackTransitions[, is_last_event := seq_len(.N) == .N, by = id]
events_withBackTransitions[, need_adjustment := id %in% adjust_ids]
events_withBackTransitions[is_last_event == TRUE, to := ifelse(
  need_adjustment,
  3,
  ifelse(from %in% c(0, 1, 2), "cens", NA)
)]
events_withBackTransitions[, transition := paste0(from, "->", to)]
events_withBackTransitions <- events_withBackTransitions[!(from == 3 & is.na(to))]
events_withBackTransitions <- merge(events_withBackTransitions, t_eskd_dt, by = 'id', all.x = TRUE)
events_withBackTransitions[need_adjustment == TRUE & is_last_event == TRUE, stop := t_eskd]
events_withBackTransitions[, episode := seq_len(.N), by = .(id, transition)]
events_withBackTransitions[, status := ifelse(to == "cens", 0, 1)]
events_withBackTransitions[, tstop2 := ifelse(
    to == "3",
    t_eskd,
    shift(tstart, type = "lead")
), by = id]
events_withBackTransitions[is.na(tstop2), tstop2 := tstop]
setorder(events_withBackTransitions, id, tstart)
events_withBackTransitions <- events_withBackTransitions[, .(id, tstart, tstop, tstop2, from, to, transition, status, episode)]

fwrite(events_withBackTransitions, file.path(dir_out, file_events_withBackTransitions_out), sep = "\t", quote = FALSE)



