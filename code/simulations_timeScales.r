library(batchtools)
library(checkmate)

# setup ----
## general
# setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
setwd("nvmetmp/wis37138/msm_kidneyFunction")
source("code/helpers_timeScales.r")
registry <- "results/registries/sim-timeScales-registry"
dir_datasets <- "/nvmetmp/wis37138/msm_kidneyFunction/results/datasets/"
repls <- 200
ncores <- 200

# simulation parameters ----

## DGP parameters and formulas ----
f_0         <- function(t) 0.02 * t^2 / (1 + 0.05 * pmax(0, t - 8)^2)
g_0         <- function(t) 0.005 * t^3
f_1         <- function(t) 0.32 * exp(-0.15 * t)
f_until_1   <- function(t) 1.30 * exp(-0.60 * t)
g_1         <- function(t) 0.14 * exp(-0.25 * t)
g_until_1   <- function(t) 0.14 * exp(-0.25 * t)

f_0_1 <- f_0
f_0_3 <- g_0
f_1_2 <- function(t) 0.48 * exp(-0.10 * t)
f_1_3 <- function(t) 0.16 * exp(-0.30 * t)

delta_0   <- -0.7 # to achieve desired type I censoring rate

beta_0_01 <- -1.8 + delta_0
beta_0_03 <- -3.9 + delta_0
beta_0_12 <- -5.5 + delta_0
beta_0_13 <- -0.6 + delta_0

beta_1_01 <- 0.0 # inspired by smoking
beta_1_03 <- 1.0
beta_1_12 <- 0.5
beta_1_13 <- 1.0

formulas_list_timeScales <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0(t) + beta_0_01 + beta_1_01 * x1
  ),
  list(
    from = 0, to = 3,
    formula = ~
      g_0(t) + beta_0_01 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(t) + f_1(t_1) + f_until_1(t_until_1) + beta_0_01 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(t) + g_1(t_1) + g_until_1(t_until_1) + beta_0_01 + beta_1_13 * x1
  )
)

formulas_list_stratified <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0_1(t) + beta_0_01 + beta_1_01 * x1
  ),
  list(
    from = 0, to = 3,
    formula = ~
      f_0_3(t) + beta_0_01 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(t) + beta_0_01 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(t) + beta_0_01 + beta_1_13 * x1
  )
)

## model formulas ----
formula_timeScales <- ped_status ~
  s(tend, by = trans_to_3) +
  s(t_1, by = trans_after_1) +
  s(t_until_1, by = trans_after_1) +
  transition * x1

formula_stratified <- ped_status ~
  s(tend, by = transition) +
  s(t_until_1, by = trans_after_1) +
  transition * x1

## other parameters ----
cut <- seq(0, 10, by = 0.05)
terminal_states <- c(2, 3)
n <- 5000
round <- 2
cens_type <- "right"
cens_dist <- "weibull"
cens_params <- c(1.5, 10.0) # shape, scale
# data <- data.frame(
#   id          = seq_len(n),
#   x1          = rbinom(n, 1, 0.5),
#   from        = 0L,
#   t           = 0 # “calendar” clock
# )
# df <- sim_pexp_msm(formulas_list_stratified, data, cut, terminal_states, round = round, add_counterfactuals = FALSE)
# table(df %>% filter(status == 1) %>% pull(transition))
# df_cens <- df %>%
#       add_censoring(type = cens_type, distribution = cens_dist, parameters = cens_params, round = round)
# table(df_cens %>% filter(status == 1) %>% pull(transition))
sim_df <- wrapper_sim(
  data = NULL, job = NULL, formulas_list_timeScales, terminal_states, cut, n = n, round = round, cens_type = cens_type,
  cens_dist = cens_dist, cens_params = cens_params)
a <- wrapper_fixedEffects(data = NULL, job = NULL, instance = sim_df, formula = formula_timeScales)

# experiment ----
if (test_directory_exists(registry)) {
  unlink(registry, recursive = TRUE)
}
if (!test_directory_exists(registry)) {
  reg <- makeExperimentRegistry(
    registry,
    packages = c("mgcv", "dplyr", "tidyr", "survival", "pammtools"),
    seed     = 11022022)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })

  addProblem(name = "sim_timeScales", fun = wrapper_sim)
  addProblem(name = "sim_stratified", fun = wrapper_sim)
  # addAlgorithm(name = "baselineHazards_timeScales", fun = wrapper_baselineHazards)
  # addAlgorithm(name = "baselineHazards_stratified", fun = wrapper_baselineHazards)
  addAlgorithm(name = "fixedEffects_timeScales", fun = wrapper_fixedEffects)
  addAlgorithm(name = "fixedEffects_stratified", fun = wrapper_fixedEffects)

  prob_df_timeScales <- data.frame(
    formulas_list = I(list(formulas_list_timeScales)),
    terminal_states = I(list(terminal_states)),
    cut = I(list(cut)),
    n = n,
    round = round,
    cens_type = cens_type,
    cens_dist = cens_dist,
    cens_params = I(list(cens_params))
  )
  prob_df_stratified <- data.frame(
    formulas_list = I(list(formulas_list_stratified)),
    terminal_states = I(list(terminal_states)),
    cut = I(list(cut)),
    n = n,
    round = round,
    cens_type = cens_type,
    cens_dist = cens_dist,
    cens_params = I(list(cens_params))
  )
  algo_df_timeScales <- data.frame(
    formula = I(list(formula_timeScales))
  )
  algo_df_stratified <- data.frame(
    formula = I(list(formula_stratified))
  )


  addExperiments(
    prob.designs = list(
      sim_timeScales = prob_df_timeScales,
      sim_stratified = prob_df_stratified
      ),
    algo.designs  = list(
      # baselineHazards_timeScales = algo_df_timeScales,
      # baselineHazards_stratified = algo_df_stratified,
      fixedEffects_timeScales = algo_df_timeScales,
      fixedEffects_stratified = algo_df_stratified
    ),
    repls = repls)

  submitJobs(ids = findNotDone())
  waitForJobs() # sometimes necessary to wait for jobs to finish and not be wrongly included in findNotDone()
}

reg     <- loadRegistry(registry, writeable = TRUE)
ids_res_msm <- findExperiments(algo.name = "msm")
ids_res_cor <- findExperiments(algo.name = "cor")
pars_msm    <- unwrap(getJobPars()) %>% as_tibble() %>% filter(algorithm == "msm")
pars_msm <- getJobPars() %>%
  as_tibble() %>%
  filter(algorithm == "msm") %>%
  transmute(
    job.id          = id,
    formula         = map_chr(algo.pars, "formula"),
    n               = map_int(prob.pars, "n"),
    round           = map_int(prob.pars, "round"),
    cens_type       = map_chr(prob.pars, "cens_type"),
    cens_dist       = map_chr(prob.pars, "cens_dist"),
    cens_params     = map(prob.pars, "cens_params")
  )

pars_cor    <- unwrap(getJobPars()) %>% as_tibble() %>% filter(algorithm == "cor")

res_msm <- reduceResultsDataTable(ids = ids_res_msm) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_msm, by = "job.id")
saveRDS(res_msm, paste0(dir_datasets, "sim-ieb-results_msm.rds"))

res_cor <- reduceResultsDataTable(ids = ids_res_cor) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_cor, by = "job.id")
saveRDS(res_cor, paste0(dir_datasets, "sim-ieb-results_cor.rds"))
