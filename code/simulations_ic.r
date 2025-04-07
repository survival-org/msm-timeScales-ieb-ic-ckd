library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)
library(batchtools)
library(checkmate)

# setup ----
# devtools::install_github("adibender/pammtools")
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
setwd("nvmetmp/wis37138/msm_kidneyFunction")
source("code/helpers_ic.r")

registry <- "results/sim-coverage-registry_baseline"
dir_dataset <- "results/dataset/"
dir_figures <- "results/figures/sim-coverage-results_baseline"
repls <- 500
ncores <- 200
formula_pexp <- "~ -3.5 + dgamma(t, 8, 2) * 6"
formula_weibull <- "~ -3.5"

# set.seed(11022022)
formula_weibull <- "~ -3.5"
# a <- wrapper_sim_pexp(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "equidistant", round = 2)
a_2 <- wrapper_sim_weibull(data = NULL, job = NULL, n = 500, formula = formula_weibull, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "beta", round = NULL)
# b <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 20, ic_point = "mid")
b_2 <- wrapper_pam(data = a_2, job = NULL, instance = a_2, bs = "ps", k = 20, ic_point = "mid")
# b_end <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 20, ic_point = "end")
# # c <- wrapper_cox(data = a, job = NULL, instance = a, ic_point = "exact")
# d <- wrapper_weibull(data = a, job = NULL, instance = a, ic_point = "end")
# # e <- wrapper_generalizedGamma(data = a, job = NULL, instance = a, ic_point = "end")


# coverage overall hazards ----
if (test_directory_exists(registry)) {
  unlink(registry, recursive = TRUE)
}
if (!test_directory_exists(registry)) {
  reg <- makeExperimentRegistry(
    registry,
    packages = c("mgcv", "dplyr", "tidyr", "pammtools", "mvtnorm", "rlang"),
    seed     = 11022022)
  reg <- loadRegistry(registry, writeable = TRUE)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })

  addProblem(name = "sim_pexp", fun = wrapper_sim_pexp)
  addProblem(name = "sim_weibull", fun = wrapper_sim_weibull)
  addAlgorithm(name = "pam", fun = wrapper_pam)
  addAlgorithm(name = "cox", fun = wrapper_cox)
  addAlgorithm(name = "weibull", fun = wrapper_weibull)
  addAlgorithm(name = "generalizedGamma", fun = wrapper_generalizedGamma)

  prob_df_pexp <- data.frame(
    formula = formula_pexp,
    ic_mechanism = c("beta", "uniform", "equidistant")
  )
  prob_df_weibull <- data.frame(
    formula = formula_weibull,
    ic_mechanism = c("beta", "uniform", "equidistant")
  )
  algo_df_pam <- data.frame(
    ic_point = c("mid", "end", "exact", "oracle", "mid_end")
  )
  algo_df_cox <- data.frame(
    ic_point = c("mid", "end", "exact", "oracle")
  )
  algo_df_weibull <- data.frame(
    ic_point = c("mid", "end", "exact", "oracle", "adjustment")
  )
  algo_df_generalizedGamma <- data.frame(
    ic_point = c("mid", "end", "exact", "oracle", "adjustment")
  )

  addExperiments(
    prob.designs = list(
      sim_pexp = prob_df_pexp,
      sim_weibull = prob_df_weibull),
    algo.designs  = list(
      pam = algo_df_pam,
      cox = algo_df_cox,
      weibull = algo_df_weibull
    ),
    repls = repls)

  submitJobs(ids = findNotDone())
  # waitForJobs()
}

reg     <- loadRegistry(registry, writeable = TRUE)
ids_res <- findExperiments()
pars    <- unwrap(getJobPars()) %>% as_tibble()

res     <- reduceResultsDataTable(ids=findDone(ids_res)) %>%
  as_tibble() %>%
  tidyr::unnest(cols = c(result)) %>%
  left_join(pars, by = "job.id")
saveRDS(res, path_dataset)

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism")

# coverage overall ----
res <- readRDS("results/datasets/sim-coverage-results_baseline.rds")
coverage <- calc_coverage(data = res, grouping_vars = grouping_vars)
coverage
View(coverage %>%
  filter(problem == "sim_weibull" & ic_mechanism=="beta" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism, problem, `coverage loghazard`)))

# RMSE ----
rmse <- calc_rmse(data = res, grouping_vars = grouping_vars)
rmse
View(rmse %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism, `RMSE loghazard`)))

# line plot ----
linePlot <- create_linePlot(data = res, grouping_vars = grouping_vars)
for(i in 1:length(linePlot)) {
  ggsave(filename = paste0(dir_figures, "/", names(linePlot)[i], ".png"), plot = linePlot[[i]], width = 10, height = 5)
}

# coverage and bias x1 ----
res <- readRDS("results/datasets/sim-coverage-results_covariate.rds")
coverage_x1 <- calc_coverage_beta(data = res, grouping_vars = grouping_vars)
coverage_x1

View(coverage_x1 %>%
  filter(problem == "sim_pexp" & ic_mechanism=="equidistant" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism, problem)))
