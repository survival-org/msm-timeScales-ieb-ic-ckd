library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)
library(patchwork)

# setup ----

setwd("nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd")
source("code/ukb/helpers_sim.r")
source("code/ukb/helpers_ieb.r")

registry <- "results/simulations/registries/sim-ieb-registry"
dir_datasets <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/"
dir_figures <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ieb/"
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

# beta_0_01 <- -3.9
# beta_0_03 <- -4
# beta_0_12 <- -3.4
# beta_0_13 <- -3.4

# beta_1_01 <- 0.2
# beta_1_03 <- 0.1
# beta_1_12 <- 0.2
# beta_1_13 <- 0.1

# beta_2_01 <- 0.2
# beta_2_03 <- 0.0
# beta_2_12 <- 0.2
# beta_2_13 <- 0.0

formulas_dgp_timeScales <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0(tend) + beta_0_01 + beta_1_01 * x1 + beta_2_01 * x2
  ),
  list(
    from = 0, to = 3,
    formula = ~
      g_0(tend) + beta_0_03 + beta_1_03 * x1 + beta_2_03 * x2
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_0(tend) + f_1(t_1) + f_until_1(t_until_1) + beta_0_12 + beta_1_12 * x1 + beta_2_12 * x2
  ),
  list(
    from = 1, to = 3,
    formula = ~
      g_0(tend) + g_1(t_1) + g_until_1(t_until_1) + beta_0_13 + beta_1_13 * x1 + beta_2_13 * x2
  )
)

formulas_dgp_stratified <- list(
  list(from = 0, to = 1,
    formula = ~
      f_0_1(tend) + beta_0_01 + beta_1_01 * x1 + beta_2_01 * x2
  ),
  list(
    from = 0, to = 3,
    formula = ~
      f_0_3(tend) + beta_0_03 + beta_1_03 * x1 + beta_2_03 * x2
  ),
  list(
    from = 1, to = 2,
    formula = ~
      f_1_2(tend) + f_until_1(t_until_1) + beta_0_12 + beta_1_12 * x1 + beta_2_12 * x2
  ),
  list(
    from = 1, to = 3,
    formula = ~
      f_1_3(tend) + g_until_1(t_until_1) + beta_0_13 + beta_1_13 * x1 + beta_2_13 * x2
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
bs <- "ps"
k <- 20

## static parameters ----
static_templates <- list(
  timeScales           = formulas_dgp_timeScales,
  stratified_fe        = formulas_dgp_stratified
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
      bs              = bs,
      k               = k
    )
    li
  }),
  names(static_templates)
)

## variable parameters ----
# dist_x1 <- c("rbinom(n, 1, 0.2)", "rbinom(n, 1, 0.5)", "rbinom(n, 1, 0.8)", "rnorm(n, 0, 1)", "rnorm(n, 0, 2)", "rnorm(n, 0, 5)")
# dist_x2 <- c("rbinom(n, 1, 0.2)", "rbinom(n, 1, 0.5)", "rbinom(n, 1, 0.8)", "rnorm(n, 0, 1)", "rnorm(n, 0, 2)", "rnorm(n, 0, 5)")

dist_x1 <- c("rnorm(n, 0, 5)")
dist_x2 <- c("rbinom(n, 1, 0.5)", "rnorm(n, 0, 1)")

# dist_x1 <- c("rbinom(n, 1, 0.5)", "rnorm(n, 0, 1)", "rnorm(n, 0, 5)")
# dist_x2 <- c("rbinom(n, 1, 0.5)", "rnorm(n, 0, 1)", "rnorm(n, 0, 5)")

prob_df_dist <- expand.grid(
  dist_x1 = dist_x1,
  dist_x2 = dist_x2,
  stringsAsFactors = FALSE
)

betas <- list(
  moderate = c(beta_0_01 = -3.9, beta_0_03 = -4.0,
               beta_0_12 = -3.4, beta_0_13 = -3.4,
               beta_1_01 = 0.2, beta_1_03 = 0.1,
               beta_1_12 = 0.2, beta_1_13 = 0.1,
               beta_2_01 = 0.2, beta_2_03 = 0.0,
               beta_2_12 = 0.2, beta_2_13 = 0.0),
  strong   = c(beta_0_01 = -3.9, beta_0_03 = -4.0,
               beta_0_12 = -3.4, beta_0_13 = -3.4,
               beta_1_01 = 0.5, beta_1_03 = 0.0,
               beta_1_12 = 0.5, beta_1_13 = 0.0,
               beta_2_01 = 0.5, beta_2_03 = 0.0,
               beta_2_12 = 0.5, beta_2_13 = 0.0)
)

## model formulas ----
formula_mod_timeScales <- ped_status ~
  s(tend, by = trans_to_3, bs = bs, k = k) +
  s(t_1, by = trans_after_1, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1 + transition * x2

formula_mod_timeScales_ieb <- ped_status ~
  s(tend, by = trans_to_3, bs = bs, k = k) +
  s(t_1, by = trans_after_1, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1

formula_mod_stratified <- ped_status ~
  s(tend, by = transition, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1 + transition * x2

formula_mod_stratified_ieb <- ped_status ~
  s(tend, by = transition, bs = bs, k = k) +
  s(t_until_1, by = trans_after_1, bs = bs, k = k) +
  transition * x1

formula_mod_list <- list(
  timeScales = formula_mod_timeScales,
  timeScales_ieb = formula_mod_timeScales_ieb,
  stratified = formula_mod_stratified,
  stratified_ieb = formula_mod_stratified_ieb
)

# trouble shooting ----
set.seed(123)
n <- 5000
data <- sim_statics$timeScales
# formulas_dgp <- formulas_dgp_timeScales
formulas_dgp <- formulas_dgp_stratified
prob_df <- prob_df_dist %>%
  crossing(beta_scenario = names(betas)) %>%
  mutate(
    beta_vals = map(beta_scenario, ~ betas[[.x]])
  )
# dist_x1 <- "rbinom(n, 1, 0.5)"
dist_x1 <- "rnorm(n, 0, 5)"
dist_x2 <- "rnorm(n, 0, 1)"
beta_vals <- prob_df$beta_vals[18]

instance <- wrapper_sim(
  job = NULL, data = data, formulas_dgp = formulas_dgp_timeScales, dist_x1 = dist_x1, dist_x2 = dist_x2, beta_vals = beta_vals,
  terminal_states = terminal_states, cut = cut, n = n, round = 2, cens_type = cens_type,
  cens_dist = cens_dist, cens_params = cens_params
)

instance$beta_vals
dim(instance$events)
table(instance$events %>% filter(status == 1) %>% pull(transition))
nrow(instance$events %>% filter(status == 0))
instance$events %>% filter(from ==0) %>% pull(status) %>% table()
instance$events %>% filter(from ==1) %>% pull(status) %>% table()

transitions <- c("0->1", "0->3", "1->2", "1->3")
for(trans in transitions){
  hist_events_trans <- instance$events %>%
    filter(transition == trans, status == 1) %>%
    ggplot(aes(x = tstop)) +
    geom_histogram(bins = 50, fill = "lightblue", color = "black") +
    labs(title = paste("Histogram of event times for transition")) +
    xlab("Event time (tend)")
  ggsave(
    filename = paste0("histogram_events_", trans, ".png"),
    plot   = hist_events_trans,
    width  = 8, height = 6
  )
}

ped <- instance$ped

ped$trans_after_1 <- factor(ped$trans_after_1,
                            levels = c("none", "1->2", "1->3"),
                            ordered = TRUE)

formula <- formula_mod_timeScales
mod <- bam(as.formula(formula)
                    , data = ped
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
)
summary(mod)

formula_ieb <- formula_mod_timeScales_ieb
mod_ieb <- bam(as.formula(formula_ieb)
                    , data = ped
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
)
summary(mod_ieb)

formula_s <- formula_mod_stratified
mod_s <- bam(as.formula(formula_s)
                    , data = ped
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
)
summary(mod_s)

msm <- wrapper_msm(data = data, job = NULL, instance = instance, formula = formula)
msm

msm_ieb <- wrapper_msm(data = data, job = NULL, instance = instance, formula = formula_ieb)
msm_ieb

msm_s <- wrapper_msm(data = data, job = NULL, instance = instance, formula = formula_s)
msm_s

cor <- wrapper_cor(data = data, job = NULL, instance = instance)
cor

# experiment ----
# if (test_directory_exists(registry)) {
#   unlink(registry, recursive = TRUE)
# }
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

  addProblem(
    name = "sim_timeScales",
    data = sim_statics$timeScales,
    fun  = function(job, data, dist_x1, dist_x2, beta_scenario, beta_vals) {
      wrapper_sim(
        job, data,
        dist_x1          = dist_x1,
        dist_x2          = dist_x2,
        beta_vals        = beta_vals,
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
    name = "sim_stratified",
    data = sim_statics$stratified,
    fun  = function(job, data, dist_x1, dist_x2, beta_scenario, beta_vals) {
      wrapper_sim(
        job, data,
        dist_x1          = dist_x1,
        dist_x2          = dist_x2,
        beta_vals        = beta_vals,
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

  prob_df <- prob_df_dist %>%
    crossing(beta_scenario = names(betas)) %>%
    mutate(
      beta_vals = map(beta_scenario, ~ betas[[.x]])
    )

  addAlgorithm(name = "msm", fun = wrapper_msm)
  addAlgorithm(name = "cor", fun = wrapper_cor)

  algo_df_msm <- data.frame(
    formula = vapply(
      formula_mod_list,
      function(f) paste0(deparse(f, width.cutoff = 500L), collapse = ""),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
  algo_df_cor <- data.frame()

  addExperiments(
    prob.designs = list(
      sim_timeScales = prob_df,
      sim_stratified = prob_df),
    algo.designs  = list(
      msm = algo_df_msm,
      cor = algo_df_cor
    ),
    repls = repls)

  start_time <- Sys.time()
  submitJobs(ids = findNotDone())
  waitForJobs() # sometimes necessary to wait for jobs to finish and not be wrongly included in findNotDone()
  end_time <- Sys.time()
  message("Total time taken for job submission and completion: ", round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
}

# adjusting registry ----
ids_all <- findJobs(reg = reg)
pars <- as.data.table( getJobPars(ids = ids_all, reg = reg) )
badJobs <- pars[
  # for each row, look into prob.pars[[i]] and test your three conditions
  sapply(prob.pars, function(pr) {
    pr$dist_x1       == "rnorm(n, 0, 5)" &&
    pr$dist_x2       != "rnorm(n, 0, 5)" &&
    pr$beta_scenario == "strong"
  }),
  job.id
]
removeExperiments(ids = badJobs, reg = reg)

# collect results ----
reg     <- loadRegistry(registry, writeable = TRUE)
ids_res_msm <- findExperiments(algo.name = "msm")
ids_res_cor <- findExperiments(algo.name = "cor")
pars_msm    <- unwrap(getJobPars(ids = ids_res_msm)) %>% as_tibble()
pars_cor    <- unwrap(getJobPars(ids = ids_res_cor)) %>% as_tibble()

res_msm <- reduceResultsDataTable(ids = ids_res_msm) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_msm, by = "job.id") %>%
  left_join(
    algo_df_msm %>% rownames_to_column(var = "formula_name"),
    by = c("formula" = "formula")
  ) %>%
  mutate(
    dist_x1_name = case_when(
      startsWith(dist_x1, "rbinom") ~ paste0("bernoulli", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x1) %>% trimws()),
      startsWith(dist_x1, "rnorm")  ~ paste0("normal", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x1) %>% trimws())
    ),
    dist_x2_name = case_when(
      startsWith(dist_x2, "rbinom") ~ paste0("bernoulli", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x2) %>% trimws()),
      startsWith(dist_x2, "rnorm")  ~ paste0("normal", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x2) %>% trimws())
    ),
    dist_x1_x2 = paste0(dist_x1_name, "_", dist_x2_name)
  )
saveRDS(res_msm, paste0(dir_datasets, "sim-ieb-results_msm.rds"))
res_msm <- readRDS(paste0(dir_datasets, "sim-ieb-results_msm.rds"))

res_cor <- reduceResultsDataTable(ids = ids_res_cor) %>%
  as_tibble() %>%
  unnest(cols = c(result)) %>%
  left_join(pars_cor, by = "job.id") %>%
  mutate(
    dist_x1_name = case_when(
      startsWith(dist_x1, "rbinom") ~ paste0("bernoulli", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x1) %>% trimws()),
      startsWith(dist_x1, "rnorm")  ~ paste0("normal", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x1) %>% trimws())
    ),
    dist_x2_name = case_when(
      startsWith(dist_x2, "rbinom") ~ paste0("bernoulli", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x2) %>% trimws()),
      startsWith(dist_x2, "rnorm")  ~ paste0("normal", gsub(".*,(\\s*[0-9.]+)\\)", "\\1", dist_x2) %>% trimws())
    ),
    dist_x1_x2 = paste0(dist_x1_name, "_", dist_x2_name)
  )
saveRDS(res_cor, paste0(dir_datasets, "sim-ieb-results_cor.rds"))
res_cor <- readRDS(paste0(dir_datasets, "sim-ieb-results_cor.rds"))

# summarize results ----
msm_summary <- res_msm %>%
  group_by(transition, problem, formula_name, dist_x1_x2) %>%
  summarise(
    coef_mean = mean(coefficient),
    bias_mean = mean(bias),
    coverage_mean = mean(coverage)
  )

cor_summary <- res_cor %>%
  group_by(problem, dist_x1_x2, from) %>%
  summarise(
    rho_mean = mean(rho),
    rho_ci = list(round(quantile(rho, probs = c(0.025, 0.975)), 4))
    ) %>%
    mutate(
      rho_lower = map_dbl(rho_ci, ~ .[1]),
      rho_upper = map_dbl(rho_ci, ~ .[2]),
      rho_ci = map_chr(rho_ci, ~ paste0("(", .[1], ", ", .[2], ")")))

# plot results ----
problem_labs <- c(sim_stratified = "SSTS DGP", sim_timeScales = "MTS DGP")

distinct_dist_x1_x2 <- unique(res_cor$dist_x1_x2)

## correlation boxplots ----

df <- res_cor
df %>% pull(rho) %>% summary()
scns <- unique(df$beta_scenario)

for (dist in distinct_dist_x1_x2) {
  for (scn in scns) {
    df_filtered <- df %>%
      filter(dist_x1_x2    == dist,
             beta_scenario == scn)

    p_cor <- plot_cor(df_filtered)

    ggsave(
      filename = file.path(
        dir_figures,
        "correlations",
        paste0(
          "correlation_boxplots_",
          dist, "_", scn,
          ".png"
        )
      ),
      plot   = p_cor,
      width  = 9,
      height = 5
    )
  }
}

## coefficient and coverage boxplots ----

df   <- res_msm
df %>% filter(transition != "progression_int") %>% pull(coefficient) %>% summary()
df %>% filter(transition == "progression_int") %>% pull(coefficient) %>% summary()
scns <- unique(df$beta_scenario)

for(dist in distinct_dist_x1_x2) {
  for(scn in scns) {
    # subset by both distance‐spec and scenario
    df_sub <- df %>%
      filter(dist_x1_x2    == dist,
             beta_scenario == scn)

    # 1) coefficient boxplot
    p1 <- plot_coefs(df_sub)
    ggsave(
      filename = file.path(
        dir_figures,
        "coefficients",
        paste0("coefficient_boxplots_",
               dist, "_", scn, ".png")
      ),
      plot   = p1,
      width  = 8, height = 8
    )

    # 2) coverage barplot
    p2 <- plot_coverage(df_sub)
    ggsave(
      filename = file.path(
        dir_figures,
        "coverages",
        paste0("coverage_boxplots_",
               dist, "_", scn, ".png")
      ),
      plot   = p2,
      width  = 8, height = 8
    )
  }
}

## facet plots ----

### correlation facet boxplot ----

# 1. Filter data for the LEFT plot
df_left <- res_cor %>%
  filter(dist_x1_x2    == "bernoulli0.5_bernoulli0.5",
         beta_scenario == "moderate")

# 2. Create the LEFT plot and add the plotmath title
p_left <- plot_cor(df_left, font_size = 20) +
  labs(title = expression(atop("Moderate Effect Size",
                               x[1] ~ " ~ Ber(0.5) & " ~ x[2] ~ " ~ Ber(0.5)")))

# 3. Filter data for the RIGHT plot
df_right <- res_cor %>%
  filter(dist_x1_x2    == "normal5_normal5",
         beta_scenario == "strong")

# 4. Create the RIGHT plot and add the plotmath title
p_right <- plot_cor(df_right, font_size = 20) +
  labs(title = expression(atop("Strong Effect Size",
                               x[1] ~ " ~ N(0,5) & " ~ x[2] ~ " ~ N(0,5)"))) +
  # Remove y-axis title, text, and ticks from the right plot
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

# 5. Combine the two plots using patchwork
combined_plot <-
  # Add a 0.5 cm margin to the RIGHT of the left plot
  (p_left + theme(plot.margin = margin(r = .5, unit = "cm"))) +

  # Add a 0.5 cm margin to the LEFT of the right plot
  (p_right + theme(plot.margin = margin(l = .5, unit = "cm"))) +

  # Add patchwork layout options
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(dir_figures, "correlations", "ieb_cor.png"),
  plot   = combined_plot, width = 12, height = 6, dpi = 300
)

### coefficient facet boxplot ----

# 1. Get the LISTS of plots from your function
# This returns a list of 2 ggplot objects for each call
df_left <- res_msm %>%
  filter(dist_x1_x2    == "bernoulli0.5_bernoulli0.5",
         beta_scenario == "moderate")

df_right <- res_msm %>%
  filter(dist_x1_x2    == "normal5_normal5",
         beta_scenario == "strong")

p_list_left <- plot_coefs(df_left, font_size = 20)
p_list_right <- plot_coefs(df_right, font_size = 20)

# 2. Add overall titles as separate "plots"
title_left_text <- expression(atop("Moderate Effect Size", x[1] ~ " ~ Ber(0.5) & " ~ x[2] ~ " ~ Ber(0.5)"))
title_left <- ggplot() + labs(title = title_left_text) + theme_void() +
              theme(plot.title = element_text(hjust = 0.5, size = 22))

title_right_text <- expression(atop("Strong Effect Size", x[1] ~ " ~ N(0,5) & " ~ x[2] ~ " ~ N(0,5)"))
title_right <- ggplot() + labs(title = title_right_text) + theme_void() +
               theme(plot.title = element_text(hjust = 0.5, size = 22))

# 3. Apply theme modifications to the individual right-side plots
p_list_right <- lapply(p_list_right, function(p) {
  p + theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
})

# 4. Assemble the final layout from all the individual pieces
# This layout explicitly defines every component's position.
layout <- "
  A#B
  C#D
  E#F
"

combined_plot <- title_left + title_right +
                 p_list_left[[1]] + p_list_right[[1]] +
                 p_list_left[[2]] + p_list_right[[2]] +
  plot_layout(
    design = layout,
    guides = "collect",
    # Heights corresponds to the 3 rows in the design
    heights = unit(c(0.1, 5, 5), c('cm', 'null', 'null')),
    # Widths now corresponds to the 3 columns: [Plots | Spacer | Plots]
    # To move panels closer, make the middle number smaller.
    widths = unit(c(1, 0.1, 1), c('null', 'cm', 'null'))
  ) &
  theme(legend.position = "bottom")

png(file.path(dir_figures, "coefficients", "ieb_coef.png"), width = 18, height = 10, units = "in", res = 300)
print(combined_plot)
dev.off()

# tables ----

dir_tables <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/tables/ieb/"

## correlation ----
processed_df_cor <- process_correlation_summary(cor_summary)
latex_cor <- convert_to_latex_cor(processed_df_cor, colsep = "2pt")
writeLines(latex_cor, file.path(dir_tables, "ieb_cor.tex"))

## bias ----
ieb_coef_bias_cov_summary <- create_csv_summary(res_msm)
View(ieb_coef_bias_cov_summary)

write.csv(ieb_coef_bias_cov_summary, file.path(dir_tables, "ieb_coef_bias_cov_summary.csv"), row.names = FALSE)
