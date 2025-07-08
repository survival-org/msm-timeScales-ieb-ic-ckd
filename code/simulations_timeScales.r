library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)

# setup ----

# setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
setwd("nvmetmp/wis37138/msm_kidneyFunction")
source("code/helpers_timeScales.r")
registry <- "results/simulations/registries/sim-timeScales-registry"
dir_datasets <- "/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/datasets/"
dir_figures <- "/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/figures/timeScales"
repls <- 500
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

formulas_dgp_timeScales_fe <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0(tend) + beta_0_01 + beta_1_01 * x1
  ),
  list(
    from = 0, to = 3,
    formula = ~
      g_0(tend) + beta_0_01 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(tend) + f_1(t_1) + f_until_1(t_until_1) + beta_0_01 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(tend) + g_1(t_1) + g_until_1(t_until_1) + beta_0_01 + beta_1_13 * x1
  )
)

formulas_dgp_stratified_fe <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0_1(tend) + beta_0_01 + beta_1_01 * x1
  ),
  list(
    from = 0, to = 3,
    formula = ~
      f_0_3(tend) + beta_0_01 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(tend) + beta_0_01 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(tend) + beta_0_01 + beta_1_13 * x1
  )
)

formulas_dgp_timeScales_bh <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0(tend) + beta_0_01
  ),
  list(
    from = 0, to = 3,
    formula = ~
      g_0(tend) + beta_0_01
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(tend) + f_1(t_1) + f_until_1(t_until_1) + beta_0_01
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(tend) + g_1(t_1) + g_until_1(t_until_1) + beta_0_01
  )
)

formulas_dgp_stratified_bh <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0_1(tend) + beta_0_01
  ),
  list(
    from = 0, to = 3,
    formula = ~
      f_0_3(tend) + beta_0_01
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(tend) + beta_0_01
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(tend) + beta_0_01
  )
)

## other parameters ----
cut <- seq(0, 10, by = 0.1)
terminal_states <- c(2, 3)
n <- 5000
round <- 1
cens_type <- "right"
cens_dist <- "weibull"
cens_params <- c(1.5, 10.0) # shape, scale
ci <- TRUE

## model formulas ----
formula_mod_timeScales_fe <- ped_status ~
  s(tend, by = trans_to_3) +
  s(t_1, by = trans_after_1) +
  s(t_until_1, by = trans_after_1) +
  transition * x1

formula_mod_stratified_fe <- ped_status ~
  s(tend, by = transition) +
  s(t_until_1, by = trans_after_1) +
  transition * x1

formula_mod_timeScales_bh <- ped_status ~
  s(tend, by = trans_to_3) +
  s(t_1, by = trans_after_1) +
  s(t_until_1, by = trans_after_1) +
  transition

formula_mod_stratified_bh <- ped_status ~
  s(tend, by = transition) +
  s(t_until_1, by = trans_after_1) +
  transition

## static parameters ----
static_templates <- list(
  timeScales_fe           = formulas_dgp_timeScales_fe,
  stratified_fe           = formulas_dgp_stratified_fe,
  timeScales_bh           = formulas_dgp_timeScales_bh,
  stratified_bh           = formulas_dgp_stratified_bh
)

sim_statics <- setNames(
  lapply(names(static_templates), function(nm) {
    fmt <- static_templates[[nm]]
    li  <- list(
      formulas_dgp    = fmt,
      terminal_states = terminal_states,
      cut             = cut,
      n               = n,
      round           = round,
      cens_type       = cens_type,
      cens_dist       = cens_dist,
      cens_params     = cens_params
    )
    if (grepl("_bh$", nm)) {
      li$ci <- ci
    }
    li
  }),
  names(static_templates)
)

# data <- data.frame(
#   id          = seq_len(n),
#   x1          = rbinom(n, 1, 0.5),
#   from        = 0L,
#   t           = 0 # “calendar” clock
# )

# df <- sim_pexp_msm(formulas_dgp_stratified, data, cut, terminal_states, round = round, add_counterfactuals = FALSE)
# table(df %>% filter(status == 1) %>% pull(transition))
# df_cens <- df %>%
#       add_censoring(type = cens_type, distribution = cens_dist, parameters = cens_params, round = round)
# table(df_cens %>% filter(status == 1) %>% pull(transition))
# sim_df <- wrapper_sim(
#   data = sim_statics$timeScales_bh, job = NULL, formulas_dgp_timeScales_bh, terminal_states, cut, n = n, round = round, cens_type = cens_type,
#   cens_dist = cens_dist, cens_params = cens_params)
# a <- wrapper_fe(data = sim_statics$timeScales_bh, job = NULL, instance = sim_df, formula = formula_mod_timeScales_fe)
# b <- wrapper_bh(
#   data = sim_statics$timeScales_bh, job = NULL, instance = sim_df, formula = formula_mod_timeScales_bh)

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

  ## fixed effects ----
  addProblem(
    name = "sim_timeScales_fe",
    data = sim_statics$timeScales_fe,
    fun  = function(job, data, n, round, cens_type, cens_dist) {
      wrapper_sim(
        formulas_dgp     = data$formulas_dgp,
        terminal_states  = data$terminal_states,
        cut              = data$cut,
        n                = data$n,
        round            = data$round,
        cens_type        = data$cens_type,
        cens_dist        = data$cens_dist,
        cens_params      = data$cens_params
      )
    }
  )
  addProblem(
    name = "sim_stratified_fe",
    data = sim_statics$stratified_fe,
    fun  = function(job, data, n, round, cens_type, cens_dist) {
      wrapper_sim(
        formulas_dgp     = data$formulas_dgp,
        terminal_states  = data$terminal_states,
        cut              = data$cut,
        n                = data$n,
        round            = data$round,
        cens_type        = data$cens_type,
        cens_dist        = data$cens_dist,
        cens_params      = data$cens_params
      )
    }
  )

  batchExport(
    export = list(
      formula_mod_timeScales_fe = formula_mod_timeScales_fe,
      formula_mod_stratified_fe = formula_mod_stratified_fe
    )
  )
  fe.wrapper <- function(job, data, instance, algo.name, ...) {
    checkmate::assert_choice(algo.name, c("fe_timeScales", "fe_stratified"))
    f <- switch(
      algo.name,
      fe_timeScales = formula_mod_timeScales_fe,
      fe_stratified  = formula_mod_stratified_fe
    )
    wrapper_fe(
      job      = job,
      data     = data,
      instance = instance,
      formula  = f,
      ...
    )
  }

  addAlgorithm(
    name = "fe",
    fun  = fe.wrapper
  )

  addExperiments(
    prob.designs = list(
      sim_timeScales_fe  = data.frame(),
      sim_stratified_fe  = data.frame()
    ),
    algo.designs = list(
      fe = data.frame(
        algo.name = c("fe_timeScales", "fe_stratified")
      )
    ),
    repls = repls
  )

  start_time <- Sys.time()
  submitJobs(findNotDone())
  waitForJobs()
  end_time <- Sys.time()
  print(paste("Total time for fixed effects:", end_time - start_time))

  ## baseline hazards ----
  addProblem(
    name = "sim_timeScales_bh",
    data = sim_statics$timeScales_bh,
    fun  = function(job, data, n, round, cens_type, cens_dist) {
      wrapper_sim(
        formulas_dgp    = data$formulas_dgp,
        terminal_states  = data$terminal_states,
        cut              = data$cut,
        n                = data$n,
        round            = data$round,
        cens_type        = data$cens_type,
        cens_dist        = data$cens_dist,
        cens_params      = data$cens_params
      )
    }
  )
  addProblem(
    name = "sim_stratified_bh",
    data = sim_statics$stratified_bh,
    fun  = function(job, data, n, round, cens_type, cens_dist) {
      wrapper_sim(
        formulas_dgp    = data$formulas_dgp,
        terminal_states  = data$terminal_states,
        cut              = data$cut,
        n                = data$n,
        round            = data$round,
        cens_type        = data$cens_type,
        cens_dist        = data$cens_dist,
        cens_params      = data$cens_params
      )
    }
  )

  batchExport(
    export = list(
      formula_mod_timeScales_bh = formula_mod_timeScales_bh,
      formula_mod_stratified_bh = formula_mod_stratified_bh
    )
  )
  bh.wrapper <- function(job, data, instance, algo.name, ...) {
    checkmate::assert_choice(algo.name, c("bh_timeScales", "bh_stratified"))
    f <- switch(
      algo.name,
      bh_timeScales = formula_mod_timeScales_bh,
      bh_stratified  = formula_mod_stratified_bh
    )
    wrapper_bh(
      job      = job,
      data     = data,
      instance = instance,
      formula  = f,
      ...
    )
  }

  addAlgorithm(
    name = "bh",
    fun  = bh.wrapper
  )

  addExperiments(
    prob.designs = list(
      sim_timeScales_bh  = data.frame(),
      sim_stratified_bh  = data.frame()
    ),
    algo.designs = list(
      bh = data.frame(
        algo.name = c("bh_timeScales", "bh_stratified")
      )
    ),
    repls = repls
  )

  start_time <- Sys.time()
  submitJobs(ids = findNotDone())
  waitForJobs() # sometimes necessary to wait for jobs to finish and not be wrongly included in findNotDone()
  end_time <- Sys.time()
  print(paste("Total time for baseline hazards:", end_time - start_time))
}

# collect results ----
reg <- loadRegistry(registry, writeable = TRUE)

## fixed effects ----
ids_fe <- findExperiments(algo.name = "fe", reg = reg)

pars_fe <- unwrap(getJobPars(reg = reg)) %>%
  as_tibble()

res_fe <- reduceResultsDataTable(ids = ids_fe, reg = reg) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_fe, by = "job.id")
saveRDS(res_fe, paste0(dir_datasets, "sim-timeScales-results_fe.rds"))

## baseline hazards ----
ids_bh <- findExperiments(algo.name = "bh", reg = reg)

pars_bh <- unwrap(getJobPars(reg = reg)) %>%
  as_tibble()

res_bh <- reduceResultsDataTable(ids = ids_bh, reg = reg) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_bh, by = "job.id")
saveRDS(res_bh, paste0(dir_datasets, "sim-timeScales-results_bh.rds"))
res_bh <- readRDS(paste0(dir_datasets, "sim-timeScales-results_bh.rds"))

# summarize results ----
fe_summary <- create_fe_table(res_fe, grouping_vars = c("problem", "algo.name"))
View(fe_summary)

bh_summary <- create_bh_table(res_bh, grouping_vars = c("problem", "algo.name"), time_scales = c("tend", "t_until_1"))
print(bh_summary)

cov_by_trans <- bh_summary %>%
  group_by(transition, problem, algo.name) %>%
  summarise(
    avg_loghazard_coverage   = mean(loghazard_coverage),
    avg_hazard_coverage      = mean(hazard_coverage),
    avg_cumu_hazard_coverage = mean(cumu_hazard_coverage),
    avg_trans_prob_coverage  = mean(trans_prob_coverage),
    .groups = "drop"
  )
print(cov_by_trans)

# plot results ----
grouping_vars <- c("problem", "algo.name")
scale <- "loghazard"
font_size <- 20
alpha <- 0.8
linePlots_transitions <- c("0->1", "0->3")

plots <- imap(linePlots_transitions, function(tr, idx) {
  create_bh_linePlots(
    data          = res_bh,
    trans         = tr,
    grouping_vars = grouping_vars,
    scale         = scale,
    font_size     = font_size,
    alpha         = alpha
  ) |>
    set_names(~ paste0(tr, "_", .x))   # prepend transition to each name
}) |> flatten()                         # → one flat named list

num_cores <- min(length(linePlots), parallel::detectCores())
registerDoParallel(cores = num_cores)
if (!dir.exists(dir_figures)) {
  dir.create(dir_figures, recursive = TRUE)
}
foreach(nm = names(plots), .packages = "ggplot2") %dopar% {
  ggsave(
    filename = file.path(dir_figures, paste0(nm, "_", scale, ".png")),
    plot     = plots[[nm]],
    width    = 10, height = 8
  )
}

stopImplicitCluster()
