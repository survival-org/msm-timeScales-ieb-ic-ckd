library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)
library(batchtools)

# setup ----
# devtools::install_github("adibender/pammtools")
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
source("code/helpers_ic.r")

registry_coverage <- "results/sim-coverage-registry"
repls <- 2

# set.seed(11022022)
# a <- sim_wrapper(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "uniform", round = 2)
# b <- coverage_wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 10, ic_point = "mid_end")
# c <- coverage_wrapper_cox(data = a, job = NULL, instance = a, ic_point = "true_time")
# d <- coverage_wrapper_generalizedGamma(data = a, job = NULL, instance = a, ic_point = "end")


# coverage overall hazards ----
if (checkmate::test_directory_exists(registry_coverage)) {
  unlink(registry_coverage, recursive = TRUE)
}
if (!checkmate::test_directory_exists(registry_coverage)) {
  reg <- makeExperimentRegistry(
    registry_coverage,
    packages = c("mgcv", "dplyr", "tidyr", "pammtools", "mvtnorm", "rlang"),
    seed     = 11022022)
  reg <- loadRegistry(registry_coverage, writeable = TRUE)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 150)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })
  addProblem(name = "coverage", fun = sim_wrapper)
  addAlgorithm(name = "pam", fun = wrapper_pam)
  addAlgorithm(name = "cox", fun = wrapper_cox)
  addAlgorithm(name = "generalizedGamma", fun = wrapper_generalizedGamma)

  prob_df <- data.frame(
    formula = "~ -3.5 + dgamma(t, 8, 2) * 6 - 1.3 * x1 + sqrt(x2)",
    round = 2,
    ic_mechanism = c("beta", "uniform")
  )
  algo_df_pam <- data.frame(
    ic_point = c("mid", "end", "true_time", "oracle", "mid_end")
  )
  algo_df_cox <- data.frame(
    ic_point = c("mid", "end", "true_time", "oracle")
  )
  algo_df_generalizedGamma <- data.frame(
    ic_point = c("mid", "end", "true_time", "oracle", "adjustment")
  )

  addExperiments(
    prob.designs = list(coverage = prob_df),
    algo.designs  = list(
      pam = algo_df_pam,
      cox = algo_df_cox,
      generalizedGamma = algo_df_generalizedGamma),
    repls = repls)

  submitJobs(findNotDone())
     # waitForJobs()
}

reg     <- loadRegistry(registry_coverage, writeable = TRUE)
ids_res <- findExperiments(prob.name = "coverage")
pars    <- unwrap(getJobPars()) %>% as_tibble()
res     <- reduceResultsDataTable(ids=findDone(ids_res)) %>%
  as_tibble() %>%
  tidyr::unnest(cols = c(result)) %>%
  left_join(pars, by = "job.id")

# coverage overall ----
coverage <- calc_coverage(data = res, grouping_vars = c("ic_point", "algorithm", "ic_mechanism"))
coverage

# RMSE ----
rmse <- calc_rmse(data = res, grouping_vars = c("ic_point", "algorithm", "ic_mechanism"))
rmse

# line plot ----
linePlot <- create_linePlot(data = res, grouping_vars = c("ic_point", "algorithm", "ic_mechanism"))
linePlot

# recovery baseline hazard ----
## TBD!!!

# recovery specific effects (linear) ----
## TBD!!!

# recovery specific effects (non-linear) ----
## TBD!!!
