library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)
library(survival)
library(mgcv)
library(pammtools)
library(batchtools)

# setup ----
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction")
source("code/helpers.r")

registry_coverage <- "output/sim-coverage-registry"
repls <- 10

# coverage overall hazards ----
if (!checkmate::test_directory_exists(registry_coverage)) {
  reg <- makeExperimentRegistry(
    "C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/output/sim-conf-int-registry",
    packages = c("mgcv", "dplyr", "tidyr", "pammtools"),
    seed     = 20052018)
  reg <- loadRegistry("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/output/sim-conf-int-registry", writeable = TRUE)
  # reg$cluster.functions = makeClusterFunctionsInteractive()
  addProblem(name   = "coverage", fun = sim_wrapper)
  addAlgorithm(name = "coverage", fun = ci_wrapper)

  algo_df <- data.frame(bs = c("tp", "ps"), stringsAsFactors = FALSE)

  addExperiments(algo.design  = list(ci = algo_df), repls = repls)

  submitJobs(findNotDone())
  # waitForJobs()

}

reg     <- loadRegistry(registry_coverage, writeable = TRUE)
ids_res <- findExperiments(prob.name = "coverage", algo.name = "coverage")
pars    <- unwrap(getJobPars()) %>% as_tibble()
res     <- reduceResultsDataTable(ids=findDone(ids_res)) %>%
  as_tibble() %>%
  tidyr::unnest(cols = c(result)) %>%
  left_join(pars)

  # RMSE and coverage hazard
res %>%
  group_by(bs, method) %>%
  summarize(
    "coverage hazard" = mean(hazard),
    "coverage cumulative hazard" = mean(cumu),
    "coverage survival probability" = mean(surv)) %>%
  ungroup() %>%
  mutate_if(is.numeric, ~round(., 3)) %>%
  rename("basis" = "bs") %>%
  knitr::kable()

# pem vs. pam ----
## simulation function
sim_wrapper <- function(data, job, n = 250, time_grid = seq(0, 10, by = 0.05)) {

  # create data set with covariates
  df <- tibble::tibble(x1 = runif(n, -3, 3), x2 = runif(n, 0, 6))
  # baseline hazard
  f0 <- function(t) {dgamma(t, 8, 2) * 6}
  # define function that generates nz exposures z(t_{z,1}), ..., z(t_{z,Q})

  sim_pexp(formula = ~ -3.5 + f0(t), data = df, cut = time_grid)

}

## estimation function
pam_wrapper <- function(data, job, instance,
  cut      = NA,
  bs       = "ps",
  mod_type = c("pem", "pam") ,
  max_time = 10) {

  if(is.na(cut)) {
    cut <- NULL
  } else {
    if(cut == "rough") {
      cut <- seq(0, max_time, by = 0.5)
    } else {
      if(cut == "fine") {
        cut <- seq(0, max_time, by = 0.2)
      }
    }
  }

  ped <- as_ped(data = instance, formula = Surv(time, status) ~ ., cut = cut, id="id")

  form <- "ped_status ~ s(tend) + s(x1) + s(x2)"
  if(mod_type == "pem") {
    form     <- ped_status ~ interval
    time_var <- "interval"
  } else {
    form     <- ped_status ~ s(tend, bs = bs, k = 10)
    time_var <- "tend"
  }

  mod <- gam(formula = form, data = ped, family = poisson(), offset = offset, method = "REML")
  # summary(mod)

  make_newdata(ped, tend=unique(tend)) %>%
    add_hazard(mod, type="link", se_mult = qnorm(0.975), time_var = time_var) %>%
    mutate(truth = -3.5 + dgamma(tend, 8, 2) * 6)

}

if(!checkmate::test_directory_exists("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/output/sim-pem-vs-pam-registry")) {
  reg <- makeExperimentRegistry("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/output/sim-pem-vs-pam-registry",
    packages = c("mgcv", "dplyr", "tidyr", "pammtools"),
    seed     = 20052018)
  # reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 2)
  addProblem(name   = "pem-vs-pam", fun = sim_wrapper)
  addAlgorithm(name = "pem-vs-pam", fun = pam_wrapper)

  algo_df <- tidyr::crossing(
    cut = c(NA, "fine", "rough"),
    mod_type = c("pem", "pam"))

  addExperiments(algo.design  = list("pem-vs-pam" = algo_df), repls = 10)
  submitJobs()
  waitForJobs()
}

reg     <- loadRegistry("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/output/sim-pem-vs-pam-registry", writeable = TRUE)
ids_pam <- findExperiments(prob.name="pem-vs-pam", algo.name="pem-vs-pam")
pars    <- unwrap(getJobPars()) %>% as_tibble()
res     <- reduceResultsDataTable(ids=findDone(ids_pam)) %>%
  as_tibble() %>%
  tidyr::unnest() %>%
  left_join(pars) %>%
  mutate(cut = case_when(is.na(cut) ~ "default", TRUE ~ cut))

res %>%
  mutate(
    sq_error = (truth - hazard)^2,
    covered = (truth >= ci_lower) & (truth <= ci_upper)) %>%
  group_by(job.id, mod_type, cut) %>%
  summarize(
    RMSE = sqrt(mean(sq_error)),
    coverage = mean(covered)) %>%
  group_by(mod_type, cut) %>%
  summarize(
    RMSE     = mean(RMSE),
    coverage = mean(coverage))

ggplot(res, aes(x=tend, y = hazard)) +
  geom_step(aes(group = job.id), alpha = 0.3) +
  geom_line(aes(y = truth, col = "truth"), lwd = 2) +
  facet_grid(cut ~ mod_type) +
  coord_cartesian(ylim=c(-5, 0)) +
  geom_smooth(aes(col="average estimate"), method="gam", formula = y ~ s(x),
    se=FALSE) +
  scale_color_brewer("", palette = "Dark2") +
  xlab("time")
