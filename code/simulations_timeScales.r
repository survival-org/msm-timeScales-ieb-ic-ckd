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
dir_figures <- "/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/figures/timeScales/"
repls <- 500
ncores <- 200
seed <- 11022022

# simulation parameters ----

## DGP parameters and formulas ----
f_0         <- function(t) 0.10 * t^2 / (0.7 + 0.04 * pmax(0, t - 3)^3)
g_0         <- function(t) 0.15 * t^2 / (0.9 + 0.01 * pmax(0, t - 1)^3)
f_1         <- function(t) 0.32 * exp(-0.15 * t)
f_until_1   <- function(t) 2.50 * exp(-0.60 * t)
g_1         <- function(t) 0.14 * exp(-0.25 * t)
g_until_1   <- function(t) 0.14 * exp(-0.25 * t)

f_0_1 <- f_0
f_0_3 <- g_0
f_1_2 <- function(t) 0.48 * exp(-0.10 * t)
f_1_3 <- function(t) 0.16 * exp(-0.30 * t)

delta_0   <- -1.0 # to achieve desired type I censoring rate

beta_0_01 <- -2.9 + delta_0
beta_0_03 <- -3 + delta_0
beta_0_12 <- -2.4 + delta_0
beta_0_13 <- -2.4 + delta_0

beta_1_01 <- 0.2
beta_1_03 <- 0.1
beta_1_12 <- 0.2
beta_1_13 <- 0.1

formulas_dgp_timeScales_fe <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0(tend) + beta_0_01 + beta_1_01 * x1
  ),
  list(
    from = 0, to = 3,
    formula = ~
      g_0(tend) + beta_0_03 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(tend) + f_1(t_1) + f_until_1(t_until_1) + beta_0_12 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(tend) + g_1(t_1) + g_until_1(t_until_1) + beta_0_13 + beta_1_13 * x1
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
      f_0_3(tend) + beta_0_03 + beta_1_03 * x1
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(tend) + f_until_1(t_until_1) +  beta_0_12 + beta_1_12 * x1
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(tend) + g_until_1(t_until_1) + beta_0_13 + beta_1_13 * x1
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
      g_0(tend) + beta_0_03
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(tend) + f_1(t_1) + f_until_1(t_until_1) + beta_0_12
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(tend) + g_1(t_1) + g_until_1(t_until_1) + beta_0_13
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
      f_0_3(tend) + beta_0_03
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(tend) + f_until_1(t_until_1) + beta_0_01
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(tend) + f_until_1(t_until_1) + beta_0_01
  )
)

## other parameters ----
cut <- seq(0, 10, by = 0.1)
terminal_states <- c(2, 3)
n <- 5000
round <- 2
cens_type <- "right"
cens_dist <- "weibull"
cens_params <- c(1.5, 10.0) # shape, scale
bs_0 <- "ps"
k <- 20
ci <- TRUE

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
      cens_params     = cens_params,
      bs_0            = bs_0,
      k               = k,
      ci              = ci
    )
    li
  }),
  names(static_templates)
)

## model formulas ----
formula_mod_timeScales_fe <- ped_status ~
  s(tend, by = trans_to_3, bs = bs_0, k = k) +
  s(t_1, by = trans_after_1, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1

formula_mod_stratified_fe <- ped_status ~
  s(tend, by = transition, bs = bs_0, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1

formula_mod_timeScales_bh <- ped_status ~
  s(tend, by = trans_to_3, bs = bs_0, k = k) +
  s(t_1, by = trans_after_1, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition

formula_mod_stratified_bh <- ped_status ~
  s(tend, by = transition, bs = bs_0, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition

# df_cens <- df %>%
#       add_censoring(type = cens_type, distribution = cens_dist, parameters = cens_params, round = round)
# table(df_cens %>% filter(status == 1) %>% pull(transition))
# instance <- wrapper_sim(
#   data = sim_statics$timeScales_fe, job = NULL, formulas_dgp = formulas_dgp_timeScales_fe, terminal_states, cut, n = n, round = round, cens_type = cens_type,
#   cens_dist = cens_dist, cens_params = cens_params)
# fe_instance <- wrapper_fe(
#   data = sim_statics$timeScales_fe, job = NULL, instance = instance, formula = formula_mod_timeScales_fe)
# instance <- wrapper_sim(
#   data = sim_statics$timeScales_fe, job = NULL, formulas_dgp = formulas_dgp_timeScales_fe, terminal_states, cut, n = n, round = round, cens_type = cens_type,
#   cens_dist = cens_dist, cens_params = cens_params)
# a <- wrapper_fe(data = sim_statics$timeScales_fe, job = NULL, instance = instance, formula = formula_mod_timeScales_fe, bs_0 = bs_0, bs = "ps", k = k)
# b <- wrapper_bh(
#   data = sim_statics$timeScales_bh, job = NULL, instance = instance, formula = formula_mod_timeScales_bh, bs_0 = "ps", bs = "ps", k = k, ci = TRUE)

# experiment ----
if (test_directory_exists(registry)) {
  unlink(registry, recursive = TRUE)
}
if (!test_directory_exists(registry)) {
  reg <- makeExperimentRegistry(
    registry,
    packages = c("mgcv", "dplyr", "tidyr", "survival", "pammtools"),
    seed     = seed)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })

  ## fixed effects ----
  addProblem(
    name = "sim_timeScales_fe",
    data = sim_statics$timeScales_fe,
    fun  = function(job, data, n, round, cens_type, cens_dist, cens_params) {
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
    fun  = function(job, data, n, round, cens_type, cens_dist, cens_params) {
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

  fe.wrapper <- function(job, data, instance, model, bs, ...) {
    checkmate::assert_choice(model, c("algo_timeScales", "algo_stratified"))
    f <- switch(
      model,
      algo_timeScales = formula_mod_timeScales_fe,
      algo_stratified = formula_mod_stratified_fe,
      stop("Unknown model")
    )

    wrapper_fe(
      job      = job,
      data     = data,
      instance = instance,
      formula  = f,
      bs_0     = data$bs_0,
      bs       = bs,
      k        = data$k,
      ...
    )
  }

  addAlgorithm(
    name = "fe",
    fun  = fe.wrapper
  )

  algo_df_fe <- expand.grid(
    model = c("algo_timeScales", "algo_stratified"),
    bs    = c("ps", "fs"),
    stringsAsFactors = FALSE
  )

  addExperiments(
    prob.designs = list(
      sim_timeScales_fe  = data.frame(),
      sim_stratified_fe  = data.frame()
    ),
    algo.designs = list(fe = algo_df_fe),
    repls = repls
  )

  start_time <- Sys.time()
  submitJobs(findNotDone())
  waitForJobs()
  end_time <- Sys.time()
  print(paste("Total time for fixed effects (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))

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
    },
    seed = seed,
    cache = TRUE
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
    },
    seed = seed,
    cache = TRUE
  )

  batchExport(
    export = list(
      formula_mod_timeScales_bh = formula_mod_timeScales_bh,
      formula_mod_stratified_bh = formula_mod_stratified_bh
    )
  )

  bh.wrapper <- function(job, data, instance, model, bs, ...) {
    checkmate::assert_choice(model, c("algo_timeScales", "algo_stratified"))
    f <- switch(
      model,
      algo_timeScales = formula_mod_timeScales_bh,
      algo_stratified = formula_mod_stratified_bh,
      stop("Unknown model")
    )

    wrapper_bh(
      job      = job,
      data     = data,
      instance = instance,
      formula  = f,
      bs_0     = data$bs_0,
      bs       = bs,
      k        = data$k,
      ci       = data$ci,
      ...
    )
  }

  addAlgorithm(name = "bh", fun  = bh.wrapper)

  algo_df_bh <- expand.grid(
    model = c("algo_timeScales", "algo_stratified"),
    bs    = c("ps", "fs"),
    stringsAsFactors = FALSE
  )

  addExperiments(
    prob.designs = list(
      sim_timeScales_bh   = data.frame(),
      sim_stratified_bh   = data.frame()
    ),
    algo.designs = list(bh = algo_df_bh),
    repls = repls
  )

  start_time <- Sys.time()
  submitJobs(ids = findNotDone())
  waitForJobs() # sometimes necessary to wait for jobs to finish and not be wrongly included in findNotDone()
  end_time <- Sys.time()
  print(paste("Total time for baseline hazards (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))
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
  left_join(pars_fe, by = "job.id") %>%
  mutate(model = paste0(model, "_", bs))
saveRDS(res_fe, paste0(dir_datasets, "sim-timeScales-results_fe.rds"))
res_fe <- readRDS(file.path(dir_datasets, "sim-timeScales-results_fe.rds"))

## baseline hazards ----
ids_bh <- findExperiments(algo.name = "bh", reg = reg)
pars_bh <- unwrap(getJobPars(reg = reg)) %>%
  as_tibble()

res_bh <- reduceResultsDataTable(ids = ids_bh, reg = reg) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_bh, by = "job.id") %>%
  mutate(model = paste0(model, "_", bs))
saveRDS(res_bh, paste0(dir_datasets, "sim-timeScales-results_bh.rds"))
res_bh <- readRDS(file.path(dir_datasets, "sim-timeScales-results_bh.rds"))

# summarize results ----
fe_summary <- create_fe_table(res_fe, grouping_vars = c("problem", "model"))
View(fe_summary)

bh_summary <- create_bh_table(res_bh, grouping_vars = c("problem", "model"), time_scales = c("tend", "t_until_1"))
saveRDS(bh_summary, paste0(dir_datasets, "sim-timeScales-results_bh_summary.rds"))
bh_summary <- readRDS(file.path(dir_datasets, "sim-timeScales-results_bh_summary.rds"))
# print(bh_summary)

cov_by_trans <- bh_summary %>%
  group_by(transition, problem, model) %>%
  summarise(
    avg_loghazard_coverage   = mean(loghazard_coverage),
    avg_hazard_coverage      = mean(hazard_coverage),
    avg_cumu_hazard_coverage = mean(cumu_hazard_coverage),
    avg_trans_prob_coverage  = mean(trans_prob_coverage),
    .groups = "drop"
  )
View(cov_by_trans)

# plot results ----
grouping_vars  <- c("problem", "model")
scales         <- c("loghazard", "hazard", "cumu_hazard", "trans_prob")
font_size      <- 20
alpha          <- 0.8

## line plots for 0->1 and 0->3 transitions ----
line_trans     <- c("0->1", "0->3")

linePlots <- expand_grid(trans = line_trans, scale = scales) |>
  pmap(function(trans, scale) {
    create_bh_linePlots(
      data          = res_bh,
      trans         = trans,
      grouping_vars = grouping_vars,
      scale         = scale,
      font_size     = font_size,
      alpha         = alpha
    ) |>
      set_names(~ paste0(trans, "_", .x, "_", scale))
  }) |>
  flatten()

num_cores <- min(length(linePlots), parallel::detectCores())
registerDoParallel(cores = num_cores)
if (!dir.exists(file.path(dir_figures, "line_and_slice_plots"))) {
  dir.create(file.path(dir_figures, "line_and_slice_plots"), recursive = TRUE)
}
start_time <- Sys.time()
foreach(nm = names(linePlots), .packages = c("ggplot2", "ragg")) %dopar% {
  ggsave(
    filename = file.path(dir_figures, "line_and_slice_plots", paste0(nm, ".png")),
    plot     = linePlots[[nm]],
    width    = 10, height = 8,
    device   = ragg::agg_png
  )
}
end_time <- Sys.time()
print(paste("Total time for line plots (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))
stopImplicitCluster()

### CONTINUE HERE; BY CHECKING IF LINEPLOTS EXISTS

## slice plots for 1->2 and 1->3 transitions ----
slice_trans  <- c("1->2", "1->3")
slice_points <- c(2, 4, 6, 8)

slicePlots <- expand_grid(trans = slice_trans, scale = scales) |>
  pmap(function(trans, scale) {
    create_bh_slicePlots(
      data          = res_bh,
      trans         = trans,
      grouping_vars = grouping_vars,
      scale         = scale,
      font_size     = font_size,
      alpha         = alpha,
      time_axis     = "tend",
      slice_axis    = "t_until_1",
      slice_points  = slice_points,
      ncol_facets   = 2
    ) |>
      set_names(~ paste0(trans, "_", .x, "_", scale))
  }) |>
  flatten()

num_cores <- min(length(slicePlots), parallel::detectCores())
registerDoParallel(cores = num_cores)
if (!dir.exists(file.path(dir_figures, "line_and_slice_plots"))) {
  dir.create(file.path(dir_figures, "line_and_slice_plots"), recursive = TRUE)
}
start_time <- Sys.time()
foreach(nm = names(slicePlots), .packages = "ggplot2") %dopar% {

  filename <- file.path(dir_figures, "line_and_slice_plots", paste0(nm, ".png"))

  if (file.exists(filename = filename)) {
    return(NULL)
  } else {
    ggsave(
      filename = filename,
      plot     = slicePlots[[nm]],
      width    = 10, height = 8,
      device   = ragg::agg_png
  )
  }

}
end_time <- Sys.time()
print(paste("Total time for line plots (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))
stopImplicitCluster()

## fixed effect coverage plots ----

algo_levels  <- c("algo_timeScales", "algo_stratified")   # left → right
bs_levels    <- c("ps", "fs")                             # within each algo
problem_labs <- c(sim_timeScales_fe = "time-scales DGP",
                  sim_stratified_fe = "stratified DGP")

for (var in unique(res_fe$variable)) {
  p <- plot_one_fixedEffect(var)
  ggsave(
    filename = file.path(dir_figures, paste0("fixedEffect_", var, ".png")),
    plot     = p,
    width    = 9, height = 5
  )
}

## baseline hazard coverage plots ----

for(trans in unique(cov_by_trans$transition)) {
  p <- plot_one_bh_coverage(cov_by_trans, trans)
  ggsave(
    filename = file.path(dir_figures, paste0("bh_coverage_", trans, ".png")),
    plot     = p,
    width    = 10, height = 6
  )
}

# analyze results ----

events_rds <- file.path(registry, "events_avg.rds")

if (!file.exists(events_rds)) {
  ncores <-  parallel::detectCores() - 20
  # ncores <- 10
  cl      <- makeCluster(ncores)
  registerDoParallel(cl)

  ids_bh_prob <- findExperiments(prob.pattern = "sim_(timeScales|stratified)_bh", reg = reg)

  events_all <- foreach(
    jid        = ids_bh_prob$job.id,
    .combine   = dplyr::bind_rows,            # same shape as before
    .packages  = c("batchtools", "dplyr")    # loaded on each worker

  ) %dopar% {

    job     <- batchtools::makeJob(jid, reg = reg)     # read RDS once
    pb_name <- job$prob.name                           # character scalar
    ped     <- job$instance$ped

    ped %>%
      dplyr::filter(transition %in% c("0->1", "0->3")) %>%
      dplyr::mutate(problem = pb_name) %>%
      dplyr::group_by(problem, transition, tend) %>%
      dplyr::summarise(
        n_events = sum(ped_status == 1),
        .groups  = "drop"
      )
  }

  stopCluster(cl)   # free cores

  events_avg <- events_all %>%                       # ← identical to your
    group_by(problem, transition, tend) %>%          #   original code
    summarise(avg_events = mean(n_events), .groups = "drop")

  saveRDS(events_avg, events_rds)                    # cache for next runs
} else {
  events_avg <- readRDS(events_rds)
}

bias_col <- "loghazard_bias"
rmse_col <- "loghazard_rmse"

summ_long <- bh_summary                                                     %>%
  filter(transition %in% c("0->1", "0->3"))                                 %>%
  mutate(problem_model = paste(problem, model, sep = " / "))                %>%
  select(tend, transition, problem_model, !!bias_col, !!rmse_col)           %>%
  rename(avg_bias = !!bias_col, avg_rmse = !!rmse_col)                      %>%
  pivot_longer(cols = c(avg_bias, avg_rmse),
               names_to  = "metric",
               values_to = "value")

events_rel <- events_avg %>%
  group_by(problem, transition) %>%
  mutate(rel_freq = avg_events / sum(avg_events)) %>%
  ungroup()

max_lines   <- max(abs(summ_long$value), na.rm = TRUE)
max_relfreq <- max(events_rel$rel_freq)
scale_f     <- max_lines / max_relfreq * 0.9    # keep 10 % head-room

events_scaled <- events_rel %>%
  mutate(value_scaled = rel_freq * scale_f)

p_bias <- create_bias_rmse_plot(summ_long, "avg_bias", "Average bias",
                                events_scaled)
p_rmse <- create_bias_rmse_plot(summ_long, "avg_rmse", "Average RMSE",
                                events_scaled)

ggsave(file.path(dir_figures, "bias_rmse_events_0-1_0-3.png"),
       plot   = gridExtra::grid.arrange(p_bias, p_rmse, nrow = 2),
       width  = 10, height = 12)
