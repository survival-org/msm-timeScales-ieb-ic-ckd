
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/ukb/helpers_ukb.r")

library(etm)
library(mvna)
library(kmi)
library(mgcv)
library(scales)
library(nnet)

# setup ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_out <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb"
file_pheno <- "pheno_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
file_events <- "events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_events_ps <- "events_ps_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_ped_events <- "ped_events_ukb_603015_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"

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

events_ps <- readRDS(file.path(dir_data, file_events_ps))

# stepwise risk factor table ----

# df_counts <- ped_events
# idx_0 <- df_counts %>% filter(from == 0) %>% pull(id) %>% unique()
# idx_1 <- df_counts %>% filter(from == 1) %>% pull(id) %>% unique()
# idx_2 <- df_counts %>% filter(from == 2) %>% pull(id) %>% unique()

# length(idx_0) # #at-risk for 0->1
# df_counts %>% filter(transition == "0->1", ped_status == 1) %>% nrow() # #events for 0->1

# length(idx_1) # #at-risk for 1->2
# df_counts %>% filter(transition == "1->2", ped_status == 1) %>% nrow() # #events for 1->2

# length(idx_2) # #at-risk for 2->3
# df_counts %>% filter(transition == "2->3", ped_status == 1) %>% nrow() # #events for 2->3

## create table ----

# # Helper to format p‐values in scientific notation: “mantissa x 10^exponent”
# format_sci <- function(pval, digits = 2) {
#   sci_str <- formatC(pval, format = "e", digits = digits)
#   parts   <- strsplit(sci_str, "e")[[1]]
#   mant    <- parts[1]
#   expn    <- as.integer(parts[2])
#   paste0(mant, " x 10^", expn)
# }

# 1) Define the “always‐included” part of the formula:
base_terms <- "
  s(tend, bs = 'ps', k = 20, by = transition) +
  s(age_onset, bs = 'ps', k = 20, by = transition_after_onset_strat) +
  s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
  sex * transition
"

# 2) For each iteration, list just the additional “rs77924615_A_G*transition”‐type terms:
extra_terms_list <- list(
  # 1st iteration: only “G × transition”
  "rs77924615_A_G * transition",

  # 2nd iteration: add PGS × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition"),

  # 3rd iteration: add Diabetes × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition"),

  # 4th iteration: add s(sc_BMI, by = transition)
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, bs = 'ps', k = 20, by = transition)"),

  # 5th iteration: add s(sc_UACR, by = transition)
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, bs = 'ps', k = 20, by = transition)",
    "s(sc_UACR, bs = 'ps', k = 20, by = transition)"),

  # 6th iteration: add smoking × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, bs = 'ps', k = 20, by = transition)",
    "s(sc_UACR, bs = 'ps', k = 20, by = transition)",
    "smoking * transition"),

  # 7th iteration: add s(eGFRcrea, by = transition)
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, bs = 'ps', k = 20, by = transition)",
    "s(sc_UACR, bs = 'ps', k = 20, by = transition)",
    "smoking * transition",
    "s(eGFRcrea, bs = 'ps', k = 20, by = transition)")

)

# 3) Friendly scenario names (one per iteration)
scenario_names <- c(
  "G only",
  "G + PGS",
  "G + PGS + Diabetes",
  "G + PGS + Diabetes + BMI",
  "G + PGS + Diabetes + BMI + uACR",
  "G + PGS + Diabetes + BMI + uACR + Smoking",
  "G + PGS + Diabetes + BMI + uACR + Smoking + eGFR"
)

# 5) Loop over each scenario using foreach
num_cores <- 7
registerDoParallel(cores = num_cores)

# The 'results' object will be created by the output of the foreach loop
results <- foreach(i = seq_along(extra_terms_list), .packages = c("mgcv", "dplyr", "stringr")) %dopar% {
  extras <- extra_terms_list[[i]]

  # Build the right‐hand side: base_terms + extras
  rhs_terms <- paste0(base_terms, collapse = " ")
  rhs      <- paste( c(strsplit(base_terms, "\\+")[[1]] %>% trimws(), extras), collapse = " + " )

  # Create the full formula object
  fmla <- as.formula( paste("ped_status ~", rhs) )

  # Fit bam()
  pam_mod <- bam(
    formula  = fmla,
    data     = ped_events,
    family   = poisson(),
    offset   = offset,
    method   = "fREML",
    discrete = TRUE
  )

  # Summarize
  su <- summary(pam_mod)
  vcov <- vcov(pam_mod)

  # Extract coefficients and stats...
  # (All the calculation code remains identical)
  coef_main  <- su$p.coeff["rs77924615_A_G"]
  se_main    <- su$se["rs77924615_A_G"]
  p_main     <- su$p.pv["rs77924615_A_G"]

  coef_t12   <- su$p.coeff["transition1->2:rs77924615_A_G"]
  se_t12     <- su$se["transition1->2:rs77924615_A_G"]
  p_t12      <- su$p.pv["transition1->2:rs77924615_A_G"]

  coef_12    <- su$p.coeff["rs77924615_A_G"] + su$p.coeff["transition1->2:rs77924615_A_G"]
  se_12      <- sqrt(su$se["rs77924615_A_G"]^2 + su$se["transition1->2:rs77924615_A_G"]^2 + 2*vcov["rs77924615_A_G", "transition1->2:rs77924615_A_G"])
  p_12       <- 2 * (1 - pnorm(abs(coef_12 / se_12)))

  coef_t23   <- su$p.coeff["transition2->3:rs77924615_A_G"]
  se_t23     <- su$se["transition2->3:rs77924615_A_G"]
  p_t23      <- su$p.pv["transition2->3:rs77924615_A_G"]

  coef_23    <- su$p.coeff["rs77924615_A_G"] + su$p.coeff["transition2->3:rs77924615_A_G"]
  se_23      <- sqrt(su$se["rs77924615_A_G"]^2 + su$se["transition2->3:rs77924615_A_G"]^2 + 2*vcov["rs77924615_A_G", "transition2->3:rs77924615_A_G"])
  p_23       <- 2 * (1 - pnorm(abs(coef_23 / se_23)))

  # coef_main_r  <- round(coef_main, 2)
  # se_main_r    <- round(se_main,   2)
  # p_main_fmt   <- format_sci(p_main, digits = 2)

  # coef_t12_r   <- round(coef_t12, 2)
  # se_t12_r     <- round(se_t12,   2)
  # p_t12_fmt    <- format_sci(p_t12,  digits = 2)

  # coef_12_r   <- round(coef_12, 2)
  # se_12_r     <- round(se_12,   2)
  # p_12_fmt    <- format_sci(p_12,  digits = 2)

  # coef_t23_r   <- round(coef_t23, 2)
  # se_t23_r     <- round(se_t23,   2)
  # p_t23_fmt    <- format_sci(p_t23,  digits = 2)

  # coef_23_r   <- round(coef_23, 2)
  # se_23_r     <- round(se_23,   2)
  # p_23_fmt    <- format_sci(p_23,  digits = 2)

  # CHANGED: The 'tibble' is now the return value of the loop iteration,
  # instead of being assigned to an external object.
  tibble(
    Scenario           = scenario_names[i],

    G_Coef             = coef_main,
    G_SE               = se_main,
    G_p_value          = p_main,

    `t1→2:G_Coef`      = coef_t12,
    `t1→2:G_SE`        = se_t12,
    `t1→2:G_p_value`   = p_t12,

    `1→2_Coef`   = coef_12,
    `1→2_SE`     = se_12,
    `1→2_p_value` = p_12,

    `t2→3:G_Coef`      = coef_t23,
    `t2→3:G_SE`        = se_t23,
    `t2→3:G_p_value`   = p_t23,

    `2→3_Coef`   = coef_23,
    `2→3_SE`     = se_23,
    `2→3_p_value` = p_23
  )
}
# Don't forget to stop the cluster after the loop
stopImplicitCluster()

# 6) Bind all rows into the final 6×N table
final_table <- bind_rows(results) %>%
  # 2) Now pivot to long form
  pivot_longer(
    cols          = -Scenario,
    names_to      = c("Effect", "Metric"),
    # UPDATED REGEX: We now include '→' and handle the pure transition names (1→2 and 2→3)
    # The first group (Effect) must now capture G, t1→2:G, 1→2, t2→3:G, 2→3
    names_pattern = "([A-Za-z0-9→:]+)_(Coef|SE|p_value)",
    values_to     = "Value"
  ) %>%
  # 3) Spread “Metric” back out to columns
  pivot_wider(
    names_from  = Metric,
    values_from = Value
  ) %>%
  # 4) Ensure the final column ordering
  select(Scenario, Effect, Coef, SE, p_value)

# 7) Print
View(final_table)

# 8) Save the final table to a CSV file
write.xlsx(final_table, file.path(dir_out, "genetics_model_results_originalCutoffs.xlsx"), rowNames = FALSE)

## convert table to latex ----
final_table <- read.xlsx(file.path(dir_out, "genetics_model_results_originalCutoffs.xlsx"))
latex_risk_factors <- convert_to_latex_risk_factors(final_table)
writeLines(latex_risk_factors, file.path(dir_out, "tables", "ukb-models-risk-factors.tex"))

# risk factor distribution table ----

risk_factor_distributions <- create_summary_for_latex(events_ps)
latex_risk_factor_distributions <- convert_to_latex_risk_factor_distributions(risk_factor_distributions)
writeLines(latex_risk_factor_distributions, file.path(dir_out, "tables", "ukb-risk-factor-distributions.tex"))

# adjustment versus weighting table ----

ped_events_with_weights <- ped_events %>%
  left_join(
    events_ps %>% select(id, from_num, stabilized_weight),
    by = c("id", "from" = "from_num")
  )

formulas <- list(
  unadjusted = "ped_status ~ s(tend, bs = 'ps', k = 20, by = transition) +
  s(age_onset, bs = 'ps', k = 20, by = transition_after_onset_strat) +
  s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
  sex * transition + rs77924615_A_G * transition",
  adjusted = "ped_status ~ s(tend, bs = 'ps', k = 20, by = transition) +
  s(age_onset, bs = 'ps', k = 20, by = transition_after_onset_strat) +
  s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
  sex * transition + rs77924615_A_G * transition +
  pgs_cross_594_umod * transition + diabetes * transition +
  s(sc_BMI, bs = 'ps', k = 20, by = transition) +
  s(sc_UACR, bs = 'ps', k = 20, by = transition) +
  smoking * transition +
  s(eGFRcrea, bs = 'ps', k = 20, by = transition)"
)

weights <- list(
  unweighted = NULL,
  weighted = ped_events_with_weights$stabilized_weight
)

scenarios <- expand.grid(
  formula_name = names(formulas),
  weight_name = names(weights),
  stringsAsFactors = FALSE
) %>%
  mutate(
    formula = formulas[formula_name],
    wvec = weights[weight_name]
  )

num_cores <- 4
registerDoParallel(cores = num_cores)

results_list <- foreach(i = 1:nrow(scenarios), .packages = c("mgcv", "dplyr")) %dopar% {

  current_scenario <- scenarios[i, ]
  fmla <- as.formula(current_scenario$formula[[1]])
  w <- current_scenario$wvec[[1]]

  m <- bam(
    formula = fmla,
    data = ped_events_with_weights,
    family = poisson(),
    weights = w,
    discrete = TRUE
  )

  su <- summary(m)
  vc <- vcov(m)

  term_main <- "rs77924615_A_G"
  coef_main <- su$p.coeff[term_main]
  se_main   <- su$se[term_main]
  p_main    <- su$p.pv[term_main]

  term_int_12 <- "transition1->2:rs77924615_A_G"
  coef_12 <- coef_main + su$p.coeff[term_int_12]
  se_12   <- sqrt(se_main^2 + su$se[term_int_12]^2 + 2 * vc[term_main, term_int_12])
  p_12    <- 2 * pnorm(abs(coef_12 / se_12), lower.tail = FALSE)

  term_int_23 <- "transition2->3:rs77924615_A_G"
  coef_23 <- coef_main + su$p.coeff[term_int_23]
  se_23   <- sqrt(se_main^2 + su$se[term_int_23]^2 + 2 * vc[term_main, term_int_23])
  p_23    <- 2 * pnorm(abs(coef_23 / se_23), lower.tail = FALSE)

  tibble(
    formula = current_scenario$formula_name,
    weighting = current_scenario$weight_name,
    term = c("0->1", "1->2", "2->3"),
    coef = c(coef_main, coef_12, coef_23),
    se = c(se_main, se_12, se_23),
    p = c(p_main, p_12, p_23)
  )
}

stopImplicitCluster()

results <- bind_rows(results_list)

out <- results %>%
  mutate(
    formula = ifelse(formula == "adjusted", "yes", "no"),
    weighting = ifelse(weighting == "weighted", "yes", "no"),
    coef = round(coef, 3),
    se = round(se, 3),
    p = if_else(round(as.numeric(p), 3) == 0, "$<$0.001", sprintf("%.3f", as.numeric(p))),
    term = gsub("->", " $\\\\rightarrow$ ", term)
  ) %>%
  rename(
    Adjustment = formula,
    Weighting = weighting,
    Transition = term,
    Coefficient = coef,
    SE = se,
    `P-value` = p
  )

write.csv2(out, file = file.path(dir_out, "weighting_and_adjusting_models_results_msm.csv"), row.names = FALSE)
out <- read.csv2(file = file.path(dir_out, "weighting_and_adjusting_models_results_msm.csv"))

latex_table <- kable(out, "latex", booktabs = TRUE, escape = FALSE,
      caption = "\\captionukbmodelsriskfactors",
      label = "ukb-models-risk-factors",
      align = "cccrrr",
      position = "!ht") %>%
  kable_styling(latex_options = "striped")

# output file requires manual adjustment of caption + label
writeLines(latex_table, file.path(dir_out, "tables", "ukb-adjustment-weighting.tex"))

