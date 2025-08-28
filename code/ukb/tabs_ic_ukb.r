
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/ukb/helpers_ukb.r")

library(etm)
library(mvna)
library(kmi)
library(mgcv)
library(scales)
library(nnet)
library(kableExtra)

# setup ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_out <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb"
file_pheno <- "pheno_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
file_events <- "events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_events_ps <- "events_ps_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_ped_events <- "ped_events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_ped_events_end <- "ped_events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_end_imputed_r0_ageScale.rds"

# load and create data ----
pheno <- fread(file.path(dir_data, file_pheno))

min_age <- 35
events <- readRDS(file.path(dir_data, file_events)) %>%
  mutate(
    to = ifelse(to == "death", 4, to),
    transition = paste0(from, "->", to),
    diabetes = ifelse(diabetes == 0, 0, 1)) %>%
  filter(agestart >= min_age)

ped_events <- readRDS(file.path(dir_data, file_ped_events)) %>%
  filter(tstart >= min_age) %>%
  mutate(
    to = ifelse(to == "death", 4, to),
    transition = paste0(from, "->", to) %>% as.factor(), # for pammtools::add_surv_prob
  ) %>%
  add_transVars()

ped_events_end <- readRDS(file.path(dir_data, file_ped_events_end)) %>%
  filter(tstart >= min_age) %>%
  mutate(
    to = ifelse(to == "death", 4, to),
    transition = paste0(from, "->", to) %>% as.factor(), # for pammtools::add_surv_prob
  ) %>%
  add_transVars()

ped_baseline <- ped_events %>%
  group_by(id, from) %>%
  slice_min(order_by = date_at_exam, with_ties = FALSE) %>%  # pick the single earliest row
  ungroup() %>%
  select(id, from, sc_UACR, sc_BMI, smoking)

events_ps <- readRDS(file.path(dir_data, file_events_ps))

# fit models ----
formula_stss_full <- "ped_status ~ s(tend, bs = 'ps', k = 20, by = transition) +
  s(age_onset, bs = 'ps', k = 20, by = transition_after_onset_strat) +
  s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
  sex * transition + rs77924615_A_G * transition +
  pgs_cross_594_umod * transition + diabetes * transition +
  s(sc_BMI, bs = 'ps', k = 20, by = transition) +
  s(sc_UACR, bs = 'ps', k = 20, by = transition) +
  smoking * transition +
  s(eGFRcrea, bs = 'ps', k = 20, by = transition)"

pam_stss_full <- mgcv::bam(
  formula = as.formula(formula_stss_full),
  data = ped_events,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

pam_stss_full_end <- mgcv::bam(
  formula = as.formula(formula_stss_full),
  data = ped_events_end,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

# define risk factors for table ----
risk_factors <- c("rs77924615_A_G", "diabetes", "smoking")

# create latex table ----
models <- list(
  mid = pam_stss_full,
  end = pam_stss_full_end
)

results_long <- purrr::map_dfr(
  names(models),
  function(model_name) {
    purrr::map_dfr(
      risk_factors,
      function(risk_factor) {
        tryCatch({
          extract_risk_factor_effects(models[[model_name]], risk_factor) %>%
            mutate(risk_factor = risk_factor)
        }, error = function(e) {
          message("Error processing ", risk_factor, " in model ", model_name, ": ", e$message)
          return(NULL)
        })
      }
    ) %>%
      mutate(estimation_point = model_name)
  }
)

results_wide <- results_long %>%
  pivot_wider(
    names_from = estimation_point,
    values_from = c(coef, se, p)
  )

risk_factor_display_names <- c(
  "rs77924615_A_G" = "G",
  "diabetes"       = "Diabetes",
  "smoking"        = "Smoking"
)

results_for_table <- results_wide %>%
  arrange(factor(risk_factor, levels = names(risk_factor_display_names)), term) %>%
  mutate(
    Transition_fmt = paste0("\\hspace{1em}", str_replace(term, "->", "$\\\\rightarrow$")),
    coef_mid_fmt   = sprintf("%.3f", coef_mid),
    se_mid_fmt     = sprintf("%.3f", se_mid),
    p_mid_fmt      = if_else(round(as.numeric(p_mid), 3) == 0, "$<$0.001", sprintf("%.3f", as.numeric(p_mid))),
    coef_end_fmt   = sprintf("%.3f", coef_end),
    se_end_fmt     = sprintf("%.3f", se_end),
    p_end_fmt      = if_else(round(as.numeric(p_end), 3) == 0, "$<$0.001", sprintf("%.3f", as.numeric(p_end)))
  ) %>%
  mutate(
    is_shaded = (row_number() %% 2 == 1),
    across(
      .cols = ends_with("_fmt"),
      .fns = ~ if_else(
        is_shaded,
        paste0("\\cellcolor{gray!10}{", . , "}"),
        .
      )
    )
  ) %>%
  select(
    risk_factor,
    Transition = Transition_fmt,
    coef_mid = coef_mid_fmt,
    se_mid = se_mid_fmt,
    p_mid = p_mid_fmt,
    coef_end = coef_end_fmt,
    se_end = se_end_fmt,
    p_end = p_end_fmt
  )

column_headers <- c(
  "Transition",
  "Coef.", "SE", "P-value",
  "Coef.", "SE", "P-value"
)

header_spec <- c(
  " " = 1,
  "Mid-Point Estimation" = 3,
  "End-Point Estimation" = 3
)

pack_rows_index <- table(factor(results_for_table$risk_factor, levels = names(risk_factor_display_names)))
names(pack_rows_index) <- risk_factor_display_names

latex_table_object <- results_for_table %>%
  select(-risk_factor) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    col.names = column_headers,
    align = "lrrrrrr",
    caption = "\\captionukbic",
    label = "ukb-ic",
    escape = FALSE
  ) %>%
  add_header_above(header_spec) %>%
  pack_rows(
    index = pack_rows_index,
    latex_align = 'l',
    bold = TRUE,
    latex_gap_space = "0.3em"
  )

latex_final_string <- latex_table_object %>%
  str_replace(
    "\\\\centering",
    "\\\\centering\n\\\\fontsize{9}{11}\\\\selectfont"
  ) %>%
  str_replace_all(
    "\\\\cmidrule\\(lr\\)",
    "\\\\cmidrule(l{3pt}r{3pt})"
  ) %>%
  str_replace(
    "\\\\begin\\{table\\}",
    "\\\\begin{table}[!h]"
  )

writeLines(latex_final_string, file.path(dir_out, "tables", "ukb-ic.tex"))
