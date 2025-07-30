library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)

# setup ----

# setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
setwd("nvmetmp/wis37138/msm_kidneyFunction")
analysis_name <- "bh"

source("code/helpers_ic.r")
registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
dir_figures <- file.path("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/figures/ic/", analysis_name)

repls <- 500
ncores <- 200
if(analysis_name == "bh") {
  formula_pexp <- "~ -3.5 + dgamma(t, 8, 2) * 6"
  formula_weibull <- "~ -3.5"
} else if(analysis_name == "fe") {
  formula_pexp <- "~ -3.5 + dgamma(t, 8, 2) * 6 -1.3*x1"
  formula_weibull <- "~ -3.5 -1.3*x1"
}

# set.seed(11022022)
# a <- wrapper_sim_pexp(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "equidistant", round = 2)
# a_2 <- wrapper_sim_weibull(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "equidistant", round = NULL)
# b <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 20, ic_point = "mid")
# b_2 <- wrapper_pam(data = a_2, job = NULL, instance = a_2, bs = "ps", k = 20, ic_point = "mid")
# b_end <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 20, ic_point = "end")
# # c <- wrapper_cox(data = a, job = NULL, instance = a, ic_point = "exact")
# d <- wrapper_weibull(data = a, job = NULL, instance = a, ic_point = "end")
# # e <- wrapper_generalizedGamma(data = a, job = NULL, instance = a, ic_point = "end")

# run ----
if (test_directory_exists(registry)) {
  unlink(registry, recursive = TRUE)
}
if (!test_directory_exists(registry)) {
  reg <- makeExperimentRegistry(
    registry,
    packages = c("mgcv", "dplyr", "tidyr", "pammtools", "mvtnorm", "rlang"),
    seed     = 11022022)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })

  addProblem(name = "sim_pexp", fun = wrapper_sim_pexp)
  addProblem(name = "sim_weibull", fun = wrapper_sim_weibull)
  addProblem(name = "sim_icenReg", fun = wrapper_sim_icenReg)
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
  prob_df_icenReg <- data.frame(
    scale = 4,
    shape = 1.5,
    inspections = 5,
    inspectLength = 2
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
  algo_df_weibull <- crossing(
    ic_point = c("mid", "end", "exact", "oracle", "adjustment"),
    fct = c("survreg", "flexsurvreg")
  )
  algo_df_generalizedGamma <- data.frame(
    ic_point = c("mid", "end", "exact", "oracle", "adjustment")
  )

  addExperiments(
    prob.designs = list(
      sim_pexp = prob_df_pexp,
      sim_weibull = prob_df_weibull,
      sim_icenReg = prob_df_icenReg),
    algo.designs  = list(
      pam = algo_df_pam,
      cox = algo_df_cox,
      weibull = algo_df_weibull,
      generalizedGamma = algo_df_generalizedGamma
    ),
    repls = repls)

  start_time <- Sys.time()
  submitJobs(ids = findNotDone())
  waitForJobs()
  end_time <- Sys.time()
  print(paste("Total time (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))
}

# save ----
analysis_name <- "fe"
registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
file_dataset <- paste0("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")

reg     <- loadRegistry(registry, writeable = TRUE)
ids_res <- findExperiments()
pars    <- unwrap(getJobPars()) %>% as_tibble()

res     <- reduceResultsDataTable(ids=findDone(ids_res)) %>%
  as_tibble() %>%
  tidyr::unnest(cols = c(result)) %>%
  left_join(pars, by = "job.id")
saveRDS(res, file_dataset)

# bh ----

## load ----
analysis_name <- "bh"
registry_bh <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
reg_bh     <- loadRegistry(registry, writeable = TRUE)
file_dataset_bh <- paste0("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
res_bh <- readRDS(file_dataset)
dir_figures_bh <- file.path("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/figures/ic/", analysis_name) # only do for analysis_name = "bh"

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism", "fct")

# analyses ----
coverage_bh <- calc_coverage(data = res_bh, grouping_vars = grouping_vars)
coverage_bh
View(coverage %>%
  filter(problem == "sim_icenReg" & ic_mechanism=="beta" & algorithm == "weibull" & ic_point %in% c("mid", "mid_end", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism, `coverage loghazard`)))

rmse <- calc_rmse(data = res, grouping_vars = grouping_vars)
rmse
View(rmse %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism, `RMSE loghazard`)))

## plots ----
scales <- c("loghazard", "hazard", "cumulativehazard", "survivalfunction")
font_size <- 20

### line plots ----
start_time <- Sys.time()
for(scale in scales) {
  linePlot <- create_linePlot(data = res_bh, grouping_vars = grouping_vars, scale = scale, font_size = font_size, alpha = 0.8)

  num_cores <- min(length(linePlot), parallel::detectCores())
  registerDoParallel(cores = num_cores)
  if (!dir.exists(file.path(dir_figures_bh, "line_plots"))) {
    dir.create(file.path(dir_figures_bh, "line_plots"), recursive = TRUE)
  }
  foreach(i = 1:length(linePlot), .packages = "ggplot2") %dopar% {
    ggsave(filename = paste0(file.path(dir_figures_bh, "line_plots"), "/", names(linePlot)[i], "_", scale, ".png"),
          plot = linePlot[[i]], width = 10, height = 8)
  }
  stopImplicitCluster()
}
end_time <- Sys.time()
print(paste("Total time for plotting (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))

### coverage barplots ----
plot_df_bh <- coverage_bh %>%
  filter((is.na(fct) | fct == "flexsurvreg")) %>%
  mutate(ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism))

for(scale in scales) {

  barplots_bh_scale <- plot_coverage_bh(data = plot_df_bh, grouping_vars = c("problem", "ic_mechanism"), scale = scale)

  dir_plots <- file.path(dir_figures_bh, "coverages", scale)
  if(!dir.exists(dir_plots)) {
    dir.create(dir_plots, recursive = TRUE)
  }
  for(i in 1:length(barplots_bh_scale)) {
    plot_name <- names(barplots_bh_scale)[i]

    ggsave(filename = paste0(dir_plots, "/", plot_name, ".png"),
          plot = barplots_bh_scale[[i]], width = 10, height = 8)
  }

}

## coverage and bias x1 ----
analysis_name <- "fe"
registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
reg     <- loadRegistry(registry, writeable = TRUE)
file_dataset <- paste0("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
dir_figures <- file.path("/nvmetmp/wis37138/msm_kidneyFunction/results/simulations/figures/ic/", analysis_name)
res <- readRDS(file_dataset)

View(coverage_x1 %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism)))

df_plot <- res %>%
  filter(
    problem != "sim_icenReg",
    algorithm != "generalizedGamma",
    (is.na(fct) | fct == "flexsurvreg")
  )

summary_fe <- summarize_fe(data = df_plot, grouping_vars = grouping_vars)
summary_fe

### coefficient boxplots ----
boxplots <- plot_coef_fe(data = df_plot, grouping_vars = c("problem", "ic_mechanism"))

if(!dir.exists(file.path(dir_figures, "coefficients"))) {
  dir.create(file.path(dir_figures, "coefficients"), recursive = TRUE)
}
for(i in 1:length(boxplots)) {
  plot_name <- names(boxplots)[i]

  ggsave(filename = paste0(dir_figures, "/coefficients/", plot_name, ".png"),
         plot = boxplots[[i]], width = 10, height = 8)
}

### coverage barplots ----
barplots <- plot_coverage_fe(data = summary_fe, grouping_vars = c("problem", "ic_mechanism"))

if(!dir.exists(file.path(dir_figures, "coverages"))) {
  dir.create(file.path(dir_figures, "coverages"), recursive = TRUE)
}
for(i in 1:length(barplots)) {
  plot_name <- names(barplots)[i]

  ggsave(filename = paste0(dir_figures, "/coverages/", plot_name, ".png"),
         plot = barplots[[i]], width = 10, height = 8)
}
