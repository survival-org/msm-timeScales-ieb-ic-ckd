library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)
library(batchtools)

# setup ----
# devtools::install_github("adibender/pammtools")
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
# setwd("nvmetmp/wis37138/msm_kidneyFunction")
source("code/helpers_ic.r")

registry_coverage <- "results/sim-coverage-registry"
repls <- 2

set.seed(11022022)
a <- sim_wrapper(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "uniform", round = 2)
b <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 10, ic_point = "mid_end")
# c <- wrapper_cox(data = a, job = NULL, instance = a, ic_point = "true_time")
d <- wrapper_weibull(data = a, job = NULL, instance = a, ic_point = "end")
# e <- wrapper_generalizedGamma(data = a, job = NULL, instance = a, ic_point = "end")


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
  addAlgorithm(name = "weibull", fun = wrapper_weibull)
  # addAlgorithm(name = "generalizedGamma", fun = wrapper_generalizedGamma)

  prob_df <- data.frame(
    formula = "~ -3.5 + dgamma(t, 8, 2) * 6 - 1.3 * x1 + sqrt(x2)",
    round = 2,
    ic_mechanism = c("beta", "uniform")
  )
  algo_df_pam <- data.frame(
    # mid_end is buggy; cannot simply use tmid for spline estimation
    # and tend for everything else (ped-splitting etc.)
    # ic_point = c("mid", "end", "true_time", "oracle", "mid_end")
    ic_point = c("mid", "end", "true_time", "oracle")
  )
  algo_df_cox <- data.frame(
    ic_point = c("mid", "end", "true_time", "oracle")
  )
  algo_df_weibull <- data.frame(
    ic_point = c("mid", "end", "true_time", "oracle", "adjustment")
  )
  # algo_df_generalizedGamma <- data.frame(
  #   ic_point = c("mid", "end", "true_time", "oracle", "adjustment")
  # )

  addExperiments(
    prob.designs = list(coverage = prob_df),
    algo.designs  = list(
      pam = algo_df_pam,
      cox = algo_df_cox,
      weibull = algo_df_weibull
      # generalizedGamma = algo_df_generalizedGamma
    ),
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
saveRDS(res, file = "results/datasets/sim-coverage-results.rds")
res <- readRDS("results/datasets/sim-coverage-results_n1000_roundNULL.rds")

grouping_vars <- c("algorithm", "ic_point", "ic_mechanism")

# coverage overall ----
coverage <- calc_coverage(data = res, grouping_vars = grouping_vars)
coverage
View(coverage %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "true_time", "adjustment")) %>%
  select(-c(ic_mechanism, `coverage loghazard`)))

# RMSE ----
rmse <- calc_rmse(data = res, grouping_vars = grouping_vars)
rmse
View(rmse %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "true_time", "adjustment")) %>%
  select(-c(ic_mechanism, `RMSE loghazard`)))

# line plot ----
linePlot <- create_linePlot(data = res, grouping_vars = grouping_vars)
linePlot
for(i in 1:length(linePlot)) {
  ggsave(filename = paste0("results/figures/overall_hazard_coverage/", names(linePlot)[i], ".png"), plot = linePlot[[i]], width = 10, height = 5)
}

# coverage and bias x1 ----
coverage_x1 <- calc_coverage_beta(data = res, grouping_vars = grouping_vars)
coverage_x1

View(coverage_x1 %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "true_time", "adjustment")) %>%
  select(-c(ic_mechanism)))

# recovery baseline hazard ----
## TBD!!!

# recovery specific effects (linear) ----
## TBD!!!

# recovery specific effects (non-linear) ----
## TBD!!!
