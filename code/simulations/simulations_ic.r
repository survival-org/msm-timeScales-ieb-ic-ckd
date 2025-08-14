library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)

# setup ----

setwd("nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd")
source("code/ukb/helpers_sim.r")
source("code/ukb/helpers_ic.r")

analysis_name <- "bh"
registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
dir_figures <- file.path("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ic/", analysis_name)
repls <- 500
ncores <- 200
seed <- 11022022

# simulation parameters ----

## DGP
if(analysis_name == "bh") {
  formula_pexp <- "~ -3.5 + dgamma(t, 8, 2) * 6"
  formula_weibull <- "~ -3.5"
} else if(analysis_name == "fe") {
  formula_pexp <- "~ -3.5 + dgamma(t, 8, 2) * 6 -1.3*x1"
  formula_weibull <- "~ -3.5 -1.3*x1"
}

## other parameters ----
cut <- seq(0, 10, by = 0.1)
n <- 500
round <- NULL

## static parameters ----
static_data  <- list(
  n               = n,
  cut             = cut,
  round           = round
)

# trouble shooting ----
set.seed(11022022)

instance <- wrapper_sim_icenReg(data = NULL, job = NULL, scale = 2, shape = 1.5, inspections = 20, inspectLength = 4)
dim(instance[[1]])
pam <- wrapper_pam(data = NULL, job = NULL, instance = instance, bs = "ps", k = 20, ic_point = "mid")

p <- ggplot(pam, aes(x = tend)) +
  geom_line(aes(y = loghazard), color = "black") +
  geom_line(aes(y = loghazard_true), color = "red")

ggsave("icenreg_pam.png", plot = p, width = 10, height = 8)

instance <- wrapper_sim_pexp(data = static_data, job = NULL, formula = formula_pexp, n = 500, cut = seq(0, 10, by = 0.1), ic = TRUE, ic_mechanism = "beta", round = 2)
instance3 <- wrapper_sim_weibull(data = NULL, job = NULL, n = 500, time_grid = seq(0, 10, by = 0.05), ic = TRUE, ic_mechanism = "equidistant", round = NULL)
b <- wrapper_pam(data = static_data, job = NULL, instance = instance, bs = "ps", k = 20, ic_point = "mid")
b_2 <- wrapper_pam(data = a_2, job = NULL, instance = a_2, bs = "ps", k = 20, ic_point = "mid")
b_end <- wrapper_pam(data = a, job = NULL, instance = a, bs = "ps", k = 20, ic_point = "end")
# c <- wrapper_cox(data = a, job = NULL, instance = a, ic_point = "exact")
d <- wrapper_weibull(data = a, job = NULL, instance = a, ic_point = "end")
# e <- wrapper_generalizedGamma(data = a, job = NULL, instance = a, ic_point = "end")

a <- res_bh %>% filter(
  problem == "sim_icenReg",
  algorithm == "weibull",
  ic_point == "adjustment"
)

mean(a$loghazard_cov)

instance <- readRDS("instance2.rds")

# run ----
if (test_directory_exists(registry)) {
  unlink(registry, recursive = TRUE)
}
if (!test_directory_exists(registry)) {
  reg <- makeExperimentRegistry(
    registry,
    packages = c("mgcv", "dplyr", "tidyr", "pammtools", "mvtnorm", "rlang"),
    seed     = seed)

  tryCatch({
    reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = ncores)
  }, error = function(e) {
    reg$cluster.functions = makeClusterFunctionsInteractive()
  })

  addProblem(
    name = "sim_pexp",
    data = static_data,
    fun  = function(job, data, formula, ic_mechanism) {
      wrapper_sim_pexp(
        data      = data,
        job       = job,
        n         = data$n,
        cut       = data$cut,
        round     = data$round,
        formula   = formula,
        ic_mechanism = ic_mechanism
      )
    }
  )

  addProblem(
    name = "sim_weibull",
    data = static_data,
    fun  = function(job, data, formula, ic_mechanism) {
      wrapper_sim_weibull(
        data      = data,
        job       = job,
        n         = data$n,
        cut       = data$cut,
        round     = data$round,
        formula   = formula,
        ic_mechanism = ic_mechanism
      )
    }
  )

  prob_df_pexp <- data.frame(
    formula = formula_pexp,
    ic_mechanism = c("beta", "uniform", "equidistant")
  )
  prob_df_weibull <- data.frame(
    formula = formula_weibull,
    ic_mechanism = c("beta", "uniform", "equidistant")
  )

  if(analysis_name == "bh") {
    addProblem(
      name = "sim_icenReg",
      data = static_data,
      fun  = function(job, data) {
        wrapper_sim_icenReg(
          data      = data,
          job       = job,
          n         = data$n,
          cut       = data$cut,
          round     = data$round
        )
      }
    )
    # prob_df_icenReg <- data.frame(
    #   scale = 4,
    #   shape = 1.5,
    #   inspections = 5,
    #   inspectLength = 2
    # )
    prob_df_icenReg <- data.frame()
  }

  addAlgorithm(name = "pam", fun = wrapper_pam)
  addAlgorithm(name = "cox", fun = wrapper_cox)
  addAlgorithm(name = "weibull", fun = wrapper_weibull)
  addAlgorithm(name = "generalizedGamma", fun = wrapper_generalizedGamma)

  algo_df_pam <- data.frame(
    ic_point = c("mid", "end", "exact")
  )
  algo_df_cox <- data.frame(
    ic_point = c("mid", "end", "exact")
  )
  algo_df_weibull <- crossing(
    ic_point = c("mid", "end", "exact", "adjustment"),
    fct = c("survreg", "flexsurvreg")
  )
  algo_df_generalizedGamma <- data.frame(
    ic_point = c("mid", "end", "exact", "adjustment")
  )

  if(analysis_name == "bh") {
    prob.designs = list(
      sim_pexp = prob_df_pexp,
      sim_weibull = prob_df_weibull,
      sim_icenReg = prob_df_icenReg
    )
  } else if(analysis_name == "fe") {
    prob.designs = list(
      sim_pexp = prob_df_pexp,
      sim_weibull = prob_df_weibull
    )
  }

  addExperiments(
    prob.designs = prob.designs,
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

# collect results ----
# analysis_name <- "bh"
# registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
file_dataset <- paste0("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")

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
reg_bh     <- loadRegistry(registry_bh, writeable = TRUE)
file_dataset_bh <- paste0("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
res_bh <- readRDS(file_dataset_bh)
dir_figures_bh <- file.path("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ic/", analysis_name) # only do for analysis_name = "bh"

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism", "fct")

## analyses ----
coverage_pw_bh <-
  calc_pointwise_coverage_bh(
    data         = res_bh,
    grouping_vars= grouping_vars,
    cut_var      = "tend",
    rounding     = 3
  )
View(coverage_pw_bh)

coverage_bh <- coverage_pw_bh %>%
  group_by(across(all_of(grouping_vars))) %>%
  summarise(
    coverage_loghazard = mean(coverage_loghazard),
    coverage_hazard = mean(coverage_hazard),
    coverage_cumu = mean(coverage_cumu),
    coverage_surv = mean(coverage_surv),
    coverage_loghazard_lower = mean(coverage_loghazard_lower),
    coverage_loghazard_upper = mean(coverage_loghazard_upper),
    coverage_hazard_lower = mean(coverage_hazard_lower),
    coverage_hazard_upper = mean(coverage_hazard_upper),
    coverage_cumu_lower = mean(coverage_cumu_lower),
    coverage_cumu_upper = mean(coverage_cumu_upper),
    coverage_surv_lower = mean(coverage_surv_lower),
    coverage_surv_upper = mean(coverage_surv_upper)
  ) %>%
  ungroup()
View(coverage_bh)

rmse_pw_bh <- calc_pointwise_rmse_bh(data = res_bh, grouping_vars = grouping_vars, cut_var = "tend", rounding = 3)
View(rmse_pw_bh)

rmse_bh <- rmse_pw_bh %>%
  group_by(across(all_of(grouping_vars))) %>%
  summarise(
    rmse_loghazard = mean(rmse_loghazard),
    rmse_hazard = mean(rmse_hazard),
    rmse_cumu = mean(rmse_cumu),
    rmse_surv = mean(rmse_surv)
  ) %>%
  ungroup()

## plots ----
scales <- c("loghazard", "hazard", "cumulativehazard", "survivalfunction")
font_size <- 24

### line plots ----
start_time <- Sys.time()
for (scale in scales) {
  linePlots <- create_linePlot(data = res_bh, grouping_vars = grouping_vars,
                               scale = scale, font_size = font_size, alpha = 0.8)

  num_cores <- min(length(linePlots), parallel::detectCores())
  registerDoParallel(cores = num_cores)

  # (1) ensure base and scale-level folders exist
  base_dir  <- file.path(dir_figures_bh, "line_plots")
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)
  scale_dir <- file.path(base_dir, scale)
  if (!dir.exists(scale_dir)) dir.create(scale_dir, recursive = TRUE)

  foreach(i = seq_along(linePlots), .packages = "ggplot2") %dopar% {
    nm <- names(linePlots)[i]

    # (2) parse name → problem and algo subfolders
    parts   <- strsplit(nm, "_", fixed = TRUE)[[1]]
    problem <- if (length(parts) >= 2) paste(parts[1:2], collapse = "_") else parts[1]
    algo    <- if (length(parts) >= 3) parts[3] else "unknown_algo"

    out_dir <- file.path(dir_figures_bh, "line_plots", scale, problem, algo)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    # (3) save into scale/problem/algo/
    ggsave(
      filename = file.path(out_dir, paste0(nm, "_", scale, ".png")),
      plot = linePlots[[i]], width = 10, height = 8
    )
  }
  stopImplicitCluster()
}
end_time <- Sys.time()
print(paste("Total time for plotting (in minutes):", as.numeric(difftime(end_time, start_time, units = "mins"))))

#### comparison of selected line plots ----
scale <- "loghazard"
linePlots_loghazards <- create_linePlot(
  data = res_bh, grouping_vars = grouping_vars,
  scale = scale, font_size = font_size, alpha = 0.8,
  show_title = FALSE   # remove titles (point 2)
)

top_name    <- "sim_pexp_pam_mid_beta_NA"
bottom_name <- "sim_pexp_pam_end_beta_NA"

p_top    <- linePlots_loghazards[[top_name]]
p_bottom <- linePlots_loghazards[[bottom_name]]

# Hide x on top; keep only bottom x-axis + title "Time" (already set in function)
p_top <- p_top + theme(axis.title.x = element_blank(),
                       axis.text.x  = element_blank())

final_two <- p_top / p_bottom +
  plot_layout(guides = "collect") &
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    plot.margin       = margin(5, 5, 5, 5)
  )

ggsave(file.path(dir_figures_bh, "line_plots", scale,
                 paste0("two_panel_sim_pexp_pam_mid_and_end_beta", scale, ".png")),
       final_two, width = 10, height = 12, dpi = 300)

### coverage barplots ----
plot_df_bh <- coverage_bh %>%
  filter((is.na(fct) | fct == "flexsurvreg")) %>%
  mutate(
    ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism),
    ic_mechanism = factor(ic_mechanism, levels = c("beta", "uniform", "equidistant")),
    algorithm = factor(algorithm, levels = c("pam", "cox", "weibull", "generalizedGamma")),
    ic_point = factor(ic_point, levels = c("exact", "adjustment", "mid", "end"))
  )

for(scale in scales) {

  barplots_bh_scale <- plot_coverage_bh(data = plot_df_bh, grouping_vars = c("problem", "ic_mechanism"),
  scale = scale, font_size = font_size, dot_size = 4)

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

### coverage barplots (facet wrap) ----
library(patchwork)

for(scale in scales) {

# 1) build the per-group plots for ONE scale
plots <- plot_coverage_bh(
  data = plot_df_bh,
  grouping_vars = c("problem", "ic_mechanism"),
  scale = scale,
  font_size = font_size-2
)

# 2) helper to fetch a plot by (problem, ic_mechanism) from the named list
pick <- function(lst, problem, ic) {
  nm <- paste(problem, ic, sep = ".")   # names come from split(..., data[grouping_vars])
  p <- lst[[nm]]
  if (is.null(p)) stop("Plot not found for ", nm)
  p
}

# 3) tiny wrapper to hide x on top rows / hide legends except the last
# helper: now also hides y-axis text/title; makes titles bold size 24
tweak <- function(p, show_x = TRUE, show_legend = FALSE, title = NULL,
                  hide_y_text = FALSE, title_size = 24, y_label = NULL) {
  if (!show_x) {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank())
  }
  if (!is.null(y_label)) {
    p <- p + labs(y = y_label)
  }
  if (hide_y_text) {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y  = element_blank())
  }
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0.5))
  }
  p
}

# assemble (only calls changed: hide_y_text = TRUE on right column; titles already set)
p11 <- tweak(pick(plots, col1, rows[1]), show_x = FALSE, show_legend = FALSE,
             title = col1_title, y_label = paste0("Coverage (", rows[1], ")"))
p12 <- tweak(pick(plots, col2, rows[1]), show_x = FALSE, show_legend = FALSE,
             title = col2_title, hide_y_text = TRUE)

p21 <- tweak(pick(plots, col1, rows[2]), show_x = FALSE, show_legend = FALSE,
             y_label = paste0("Coverage (", rows[2], ")"))
p22 <- tweak(pick(plots, col2, rows[2]), show_x = FALSE, show_legend = FALSE,
             hide_y_text = TRUE)

p31 <- tweak(pick(plots, col1, rows[3]), show_x = TRUE,  show_legend = FALSE,
             y_label = paste0("Coverage (", rows[3], ")"))
p31 <- p31 + theme(axis.title.y = element_text(hjust = +0.8))
p32 <- tweak(pick(plots, col2, rows[3]), show_x = TRUE,  show_legend = TRUE,
             hide_y_text = TRUE)

final_fig <-
  (p11 | p12) /
  (p21 | p22) /
  (p31 | p32) +
  plot_layout(guides = "collect") &
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) &     # new: single legend row
  theme(
    legend.position      = "bottom",                           # changed: was coords
    legend.direction     = "horizontal",
    legend.justification = "center",
    plot.margin          = margin(5, 12, 30, 5)               # a bit more bottom space
  )

# 6) save
ggsave(file.path(dir_figures_bh, "coverages",
                 paste0("coverage_grid_", scale, ".png")),
       final_fig, width = 12, height = 12, dpi = 300)

}

# fe ----
analysis_name <- "fe"
registry <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
reg     <- loadRegistry(registry, writeable = TRUE)
file_dataset <- paste0("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
dir_figures <- file.path("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ic/", analysis_name)
res <- readRDS(file_dataset)

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism", "fct")

View(coverage_x1 %>%
  filter(ic_mechanism=="beta" & ic_point %in% c("mid", "end", "exact", "adjustment")) %>%
  select(-c(ic_mechanism)))

df_plot <- res %>%
  filter(
    problem != "sim_icenReg",
    algorithm != "generalizedGamma",
    (is.na(fct) | fct == "flexsurvreg")
  ) %>%
  mutate(
    ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism),
    ic_mechanism = factor(ic_mechanism, levels = c("beta", "uniform", "equidistant")),
    algorithm = factor(algorithm, levels = c("pam", "cox", "weibull")),
    ic_point = factor(ic_point, levels = c("exact", "adjustment", "mid", "end"))
  )

summary_fe <- summarize_fe(data = df_plot, grouping_vars = grouping_vars)
summary_fe

## coefficient boxplots ----
boxplots <- plot_coef_fe(data = df_plot, grouping_vars = c("problem", "ic_mechanism"), font_size = 24)

if(!dir.exists(file.path(dir_figures, "coefficients"))) {
  dir.create(file.path(dir_figures, "coefficients"), recursive = TRUE)
}
for(i in 1:length(boxplots)) {
  plot_name <- names(boxplots)[i]

  ggsave(filename = paste0(dir_figures, "/coefficients/", plot_name, ".png"),
         plot = boxplots[[i]], width = 10, height = 8)
}

### coefficient boxplots (facet wrap) ----

pick <- function(lst, problem, ic) {
  nm <- paste(problem, ic, sep = ".")  # Names come from split(..., data[grouping_vars])
  p <- lst[[nm]]
  if (is.null(p)) stop("Plot not found for ", nm)
  p
}

tweak <- function(p, title = NULL, ic_mechanism = NULL, show_x = TRUE, title_size = 24, y_label = NULL) {
  # Dynamically set the y-axis label with the given IC mechanism
  if (!is.null(ic_mechanism)) {
    y_label <- as.expression(bquote(hat(beta)[x[1]] ~ .(ic_mechanism)))  # Dynamic y-axis label
  }

  # Hide x-axis title and labels if show_x is FALSE
  if (!show_x) {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x  = element_blank())  # Hide x-axis labels
  }

  # Add dynamic y-axis label if provided
  if (!is.null(y_label)) {
    p <- p + labs(y = y_label)
  }

  # Add title and style it
  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(face = "bold", size = title_size, hjust = 0.5))
  }

  p
}

# Columns (DGPs):
col1 <- "sim_pexp";    col1_title <- "Piecewise Exponential DGP"
col2 <- "sim_weibull"; col2_title <- "Weibull DGP"

# Rows (IC mechanisms):
rows <- c("beta", "uniform", "equidistant")

# Row 1: (no x, no legend)
p11 <- tweak(pick(boxplots, col1, rows[1]), title = col1_title, ic_mechanism = "beta", show_x = FALSE)
p12 <- tweak(pick(boxplots, col2, rows[1]), title = col2_title, ic_mechanism = "beta", show_x = FALSE)

# Row 2: (no x, no legend)
p21 <- tweak(pick(boxplots, col1, rows[2]), title = NULL, ic_mechanism = "uniform", show_x = FALSE)
p22 <- tweak(pick(boxplots, col2, rows[2]), title = NULL, ic_mechanism = "uniform", show_x = FALSE)

# Row 3: (show x, KEEP legend on one panel)
p31 <- tweak(pick(boxplots, col1, rows[3]), title = NULL, ic_mechanism = "equidistant", show_x = TRUE)
p32 <- tweak(pick(boxplots, col2, rows[3]), title = NULL, ic_mechanism = "equidistant", show_x = TRUE)

# 5) Combine using patchwork
final_fig <-
  (p11 | p12) /
  (p21 | p22) /
  (p31 | p32) +
  plot_layout(
    guides = "collect",                  # Collects legends in a single row at the bottom
    design = "
      AA
      BB
      CC
    "
  ) &
  theme(
    legend.position      = "bottom",       # Position legend at the bottom
    legend.direction     = "horizontal",   # Legend horizontal
    legend.justification = "center",       # Centered legend
    plot.margin          = margin(5, 5, 20, 5)   # Extra space at bottom for legend
  )

# ggsave(file.path(dir_figures, "coefficients", paste0("boxplot_grid.png")),
#        final_fig, width = 12, height = 12, dpi = 300)

png(file.path(dir_figures, "coefficients", paste0("boxplot_grid.png")), width = 12, height = 12, units = "in", res = 300)
print(final_fig)
dev.off()

## coverage barplots ----
barplots <- plot_coverage_fe(data = summary_fe, grouping_vars = c("problem", "ic_mechanism"))

if(!dir.exists(file.path(dir_figures, "coverages"))) {
  dir.create(file.path(dir_figures, "coverages"), recursive = TRUE)
}
for(i in 1:length(barplots)) {
  plot_name <- names(barplots)[i]

  ggsave(filename = paste0(dir_figures, "/coverages/", plot_name, ".png"),
         plot = barplots[[i]], width = 10, height = 8)
}

# tables ----

dir_tables <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/tables/ic/"

### coverage bh table ----
plot_df_bh <- coverage_bh %>%
  filter((is.na(fct) | fct == "flexsurvreg")) %>%
  mutate(
    ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism),
    ic_mechanism = factor(ic_mechanism, levels = c("beta", "uniform", "equidistant")),
    algorithm = factor(algorithm, levels = c("pam", "cox", "weibull", "generalizedGamma")),
    ic_point = factor(ic_point, levels = c("exact", "adjustment", "mid", "end"))
  )

loghazard_summary_df <- create_coverage_summary_bh(plot_df_bh, coverage_type = "loghazard")
hazard_summary_df <- create_coverage_summary_bh(plot_df_bh, coverage_type = "hazard")
surv_summary_df <- create_coverage_summary_bh(plot_df_bh, coverage_type = "surv")
cumu_summary_df <- create_coverage_summary_bh(plot_df_bh, coverage_type = "cumu")

loghazard_latex <- generate_multi_dgp_latex(loghazard_summary_df, coverage_type = "loghazard")
hazard_latex <- generate_multi_dgp_latex(hazard_summary_df, coverage_type = "hazard")
surv_latex <- generate_multi_dgp_latex(surv_summary_df, coverage_type = "surv")
cumu_latex <- generate_multi_dgp_latex(cumu_summary_df, coverage_type = "cumu")

writeLines(loghazard_latex, file.path(dir_tables, "ic-bh-coverage-loghazard.tex"))
writeLines(hazard_latex, file.path(dir_tables, "ic-bh-coverage-hazard.tex"))
writeLines(surv_latex, file.path(dir_tables, "ic-bh-coverage-surv.tex "))
writeLines(cumu_latex, file.path(dir_tables, "ic-bh-coverage-cumu.tex"))

### bias fe table ----
bias_summary_df <- create_bias_summary_fe(summary_fe)
bias_latex <- generate_bias_latex_fe(bias_summary_df)
writeLines(bias_latex, file.path(dir_tables, "fe_bias.tex"))

### coverage fe table ----
coverage_summary_df <- create_coverage_summary_fe(summary_fe)
coverage_latex <- generate_coverage_latex_fe(coverage_summary_df)
writeLines(coverage_latex, file.path(dir_tables, "fe_coverage.tex"))
