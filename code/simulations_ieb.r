library(batchtools)
library(checkmate)

# setup ----
## general
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
# setwd("nvmetmp/wis37138/msm_kidneyFunction")
source("code/helpers_ieb.r")
registry <- "results/registries/sim-ieb-registry"
dir_datasets <- "C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/results/datasets/"
dir_figures <- "C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/results/figures/"
repls <- 500
ncores <- 200

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

  addProblem(name = "sim", fun = wrapper_sim)
  addAlgorithm(name = "msm", fun = wrapper_msm)
  addAlgorithm(name = "cor", fun = wrapper_cor)

  prob_df <- data.frame(
    cens_type = c("none")
  )
  algo_df_msm <- data.frame(
    formula = c(
      "ped_status ~ s(tend, by=transition) + transition + x1*transition",
      "ped_status ~ s(tend, by=transition) + transition + x1*transition + x2*transition"
    )
  )
  algo_df_cor <- data.frame(
  )

  addExperiments(
    prob.designs = list(sim = prob_df),
    algo.designs  = list(
      msm = algo_df_msm,
      cor = algo_df_cor
    ),
    repls = repls)

  submitJobs(ids = findNotDone())
  # waitForJobs()
}

reg     <- loadRegistry(registry, writeable = TRUE)
ids_res_msm <- findExperiments(algo.name = "msm")
ids_res_cor <- findExperiments(algo.name = "cor")
pars_msm    <- unwrap(getJobPars()) %>% as_tibble() %>% filter(algorithm == "msm")
pars_cor    <- unwrap(getJobPars()) %>% as_tibble() %>% filter(algorithm == "cor")

res_msm <- reduceResultsDataTable(ids = ids_res_msm) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_msm, by = "job.id")
saveRDS(res_msm, paste0(dir_datasets, "sim-ieb-results_msm.rds"))
saveRDS(res_msm, "/nvmetmp/wis37138/msm_kidneyFunction/results/datasets/sim-ieb-results_msm.rds")

res_cor <- reduceResultsDataTable(ids = ids_res_cor) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_cor, by = "job.id")
saveRDS(res_cor, paste0(dir_datasets, "sim-ieb-results_cor.rds"))

# interpretation ----

## coverage bias ----
res_msm <- readRDS(paste0(dir_datasets, "sim-ieb-results_msm.rds"))

coverage_bias <- calc_coverage_bias(res_msm)
coverage_bias

unadj_onset <- res_msm %>%
  filter(!grepl("x2", formula)) %>%
  filter(grepl("onset", transition))

unadj_progression <- res_msm %>%
  filter(!grepl("x2", formula)) %>%
  filter(grepl("progression", transition))

adj_onset <- res_msm %>%
  filter(grepl("x2", formula)) %>%
  filter(grepl("onset", transition))

adj_progression <- res_msm %>%
  filter(grepl("x2", formula)) %>%
  filter(grepl("progression", transition))

df_list <- list(
  unadj_onset = unadj_onset,
  unadj_progression = unadj_progression,
  adj_onset = adj_onset,
  adj_progression = adj_progression
)

for(i in 1:length(df_list)) {
  df <- df_list[[i]]
  p <- ggplot(df, aes(x = coefficient)) +
    geom_histogram(bins = 30) +
    labs(x = "Coefficient",
         y = "Frequency") +
    theme_minimal()
  ggsave(paste0(dir_figures, "sim-ieb-results/", names(df_list)[i], ".png"), p, width = 8, height = 5)
}

## correlation ----
res_cor <- readRDS(paste0(dir_datasets, "sim-ieb-results_cor.rds"))

res_cor %>% filter(var1=="x1" & var2=="x2") %>%
  group_by(from) %>%
  summarise(
    correlation = mean(rho)
  )

# testing ----
## simulation parameters
n <- 2500
formulas_list <- list(
  list(from = 0, to = 1, formula = ~ -3.5 + dgamma(t, 8, 2) * 6 + 1.3*x1 + 0.8*x2),
  list(from = 0, to = "death", formula = ~ -3.0 + dgamma(t, 8, 2) * 6 + 1.3*x1 + 0.8*x2),
  list(1, "death", ~ -2.1 + 1.3*x1 + 0.8*x2)
)
terminal_states <- c("death")
cut <- seq(0, 5, by = 0.01)
round <- 2
cens_type <- "none"

## estimation parameters
formula <- ped_status ~ s(tend, by=transition) + transition + x1*transition

set.seed(123)
a <- wrapper_sim(data = NULL, job = NULL, formulas_list = formulas_list, terminal_states = terminal_states, cut = cut, n = n, round = round, cens_type = cens_type)
