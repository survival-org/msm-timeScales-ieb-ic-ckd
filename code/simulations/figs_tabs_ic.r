library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)
library(patchwork)

# setup ----

setwd("nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd")
source("code/simulations/helpers_sim.r")
source("code/simulations/helpers_ic.r")

# bh ----

## load ----
analysis_name <- "bh"
registry_bh <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
reg_bh     <- loadRegistry(registry_bh, writeable = TRUE)
file_dataset_bh <- paste0("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
res_bh <- readRDS(file_dataset_bh)
dir_figures_bh <- file.path("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ic/", analysis_name) # only do for analysis_name = "bh"
dir_tables <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/tables/ic/"

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism", "fct")

## summarise ----
cov_pw_bh <-
  calc_pointwise_coverage_bh(
    data         = res_bh,
    grouping_vars= grouping_vars,
    cut_var      = "tend",
    rounding     = 3
  )
View(cov_pw_bh)

cov_bh <- cov_pw_bh %>%
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
View(cov_bh)

plot_df_cov_bh <- cov_bh %>%
  filter((is.na(fct) | fct == "flexsurvreg")) %>%
  mutate(
    ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism),
    ic_mechanism = factor(ic_mechanism, levels = c("beta", "uniform", "equidistant")),
    algorithm = factor(algorithm, levels = c("pam", "cox", "weibull", "generalizedGamma"), labels = c("PAM", "Cox", "Weibull AFT", "Gen. Gamma AFT")),
    ic_point = factor(ic_point, levels = c("exact", "adjustment", "mid", "end"))
  )

# rmse_pw_bh <- calc_pointwise_rmse_bh(data = res_bh, grouping_vars = grouping_vars, cut_var = "tend", rounding = 3)
# View(rmse_pw_bh)

# rmse_bh <- rmse_pw_bh %>%
#   group_by(across(all_of(grouping_vars))) %>%
#   summarise(
#     rmse_loghazard = mean(rmse_loghazard),
#     rmse_hazard = mean(rmse_hazard),
#     rmse_cumu = mean(rmse_cumu),
#     rmse_surv = mean(rmse_surv)
#   ) %>%
#   ungroup()

## plots ----
scales <- c("loghazard", "hazard", "cumulativehazard", "survivalfunction")
font_size <- 24

### line plots ----

#### all ----
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
linePlots_scale <- create_linePlot(
  data = res_bh, grouping_vars = grouping_vars,
  scale = scale, font_size = font_size, alpha = 0.8,
  show_title = FALSE   # remove titles (point 2)
)

top_name    <- "sim_pexp_pam_mid_beta_NA"
bottom_name <- "sim_pexp_pam_end_beta_NA"
p_top    <- linePlots_scale[[top_name]]
p_bottom <- linePlots_scale[[bottom_name]]
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
                 paste0("ic_bh_linePlots_", scale, ".png")),
       final_two, width = 10, height = 12, dpi = 300)

### coverage barplots ----

#### individual ----
for(scale in scales) {

  barplots_bh_scale <- plot_coverage_bh(data = plot_df_cov_bh, grouping_vars = c("problem", "ic_mechanism"),
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

#### facet wrap ----
col1 <- "sim_pexp";    col1_title <- "Piecewise Exponential DGP"
col2 <- "sim_weibull"; col2_title <- "Weibull DGP"
rows <- c("beta", "uniform", "equidistant")

for(scale in scales) {

  plots <- plot_coverage_bh(
    data = plot_df_cov_bh,
    grouping_vars = c("problem", "ic_mechanism"),
    scale = scale,
    font_size = font_size-2
  )

  p11 <- tweak_bh(pick(plots, col1, rows[1]), show_x = FALSE, show_legend = FALSE,
              title = col1_title, y_label = paste0("Coverage\n(", rows[1], ")"))
  p12 <- tweak_bh(pick(plots, col2, rows[1]), show_x = FALSE, show_legend = FALSE,
              title = col2_title, hide_y_text = TRUE)

  p21 <- tweak_bh(pick(plots, col1, rows[2]), show_x = FALSE, show_legend = FALSE,
              y_label = paste0("Coverage\n(", rows[2], ")"))
  p22 <- tweak_bh(pick(plots, col2, rows[2]), show_x = FALSE, show_legend = FALSE,
              hide_y_text = TRUE)

  p31 <- tweak_bh(pick(plots, col1, rows[3]), show_x = TRUE,  show_legend = FALSE,
              y_label = paste0("Coverage\n(", rows[3], ")"))
  p31 <- p31 + theme(axis.title.y = element_text(hjust = +0.8))
  p32 <- tweak_bh(pick(plots, col2, rows[3]), show_x = TRUE,  show_legend = TRUE,
              hide_y_text = TRUE)

  ic_bh_cov_scale <-
    (p11 | p12) /
    (p21 | p22) /
    (p31 | p32) +
    plot_layout(guides = "collect") &
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) &
    theme(
      legend.position      = "bottom",
      legend.direction     = "horizontal",
      legend.justification = "center",
      plot.margin          = margin(5, 12, 30, 5)
    )

  scale_name <- switch(
    scale,
    loghazard = "loghazard",
    hazard = "hazard",
    cumulativehazard = "cumu",
    survivalfunction = "surv"
  )
  ggsave(file.path(dir_figures_bh,
                  paste0("ic_bh_cov_", scale_name, ".png")),
        ic_bh_cov_scale, width = 12, height = 12, dpi = 300)

}

## tables ----

### coverage bh table ----
loghazard_summary_df <- create_coverage_summary_bh(plot_df_cov_bh, coverage_type = "loghazard")
hazard_summary_df <- create_coverage_summary_bh(plot_df_cov_bh, coverage_type = "hazard")
surv_summary_df <- create_coverage_summary_bh(plot_df_cov_bh, coverage_type = "surv")
cumu_summary_df <- create_coverage_summary_bh(plot_df_cov_bh, coverage_type = "cumu")

loghazard_latex <- generate_multi_dgp_latex(loghazard_summary_df, coverage_type = "loghazard")
hazard_latex <- generate_multi_dgp_latex(hazard_summary_df, coverage_type = "hazard")
surv_latex <- generate_multi_dgp_latex(surv_summary_df, coverage_type = "surv")
cumu_latex <- generate_multi_dgp_latex(cumu_summary_df, coverage_type = "cumu")

writeLines(loghazard_latex, file.path(dir_tables, "ic-bh-coverage-loghazard.tex"))
writeLines(hazard_latex, file.path(dir_tables, "ic-bh-coverage-hazard.tex"))
writeLines(surv_latex, file.path(dir_tables, "ic-bh-coverage-surv.tex "))
writeLines(cumu_latex, file.path(dir_tables, "ic-bh-coverage-cumu.tex"))

# fe ----

## load ----
analysis_name <- "fe"
registry_fe <- paste0("results/simulations/registries/sim-ic-registry_", analysis_name)
reg_fe     <- loadRegistry(registry_fe, writeable = TRUE)
file_dataset_fe <- paste0("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/sim-ic-results_", analysis_name, ".rds")
res <- readRDS(file_dataset_fe)
dir_figures_fe <- file.path("/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ic/", analysis_name)
dir_tables <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/tables/ic/"

grouping_vars <- c("problem", "algorithm", "ic_point", "ic_mechanism", "fct")

## summarise ----
plot_df_fe <- res %>%
  filter(
    problem != "sim_icenReg",
    algorithm != "generalizedGamma",
    (is.na(fct) | fct == "flexsurvreg")
  ) %>%
  mutate(
    ic_mechanism = ifelse(problem == "sim_icenReg", "uniform", ic_mechanism),
    ic_mechanism = factor(ic_mechanism, levels = c("beta", "uniform", "equidistant")),
    algorithm = factor(algorithm, levels = c("pam", "cox", "weibull"), labels = c("PAM", "Cox", "Weibull AFT")),
    ic_point = factor(ic_point, levels = c("exact", "adjustment", "mid", "end"))
  )

summary_fe <- summarize_fe(data = plot_df_fe, grouping_vars = grouping_vars)
summary_fe

## plots ----
font_size <- 24

### coefficient boxplots ----
boxplots <- plot_coef_fe(data = plot_df_fe, grouping_vars = c("problem", "ic_mechanism"), font_size = font_size)

#### individual ----
if(!dir.exists(file.path(dir_figures_fe, "coefficients"))) {
  dir.create(file.path(dir_figures_fe, "coefficients"), recursive = TRUE)
}
for(i in 1:length(boxplots)) {
  plot_name <- names(boxplots)[i]

  ggsave(filename = paste0(dir_figures_fe, "/coefficients/", plot_name, ".png"),
         plot = boxplots[[i]], width = 10, height = 8)
}

#### facet wrap ----
col1 <- "sim_pexp";    col1_title <- "Piecewise Exponential DGP"
col2 <- "sim_weibull"; col2_title <- "Weibull DGP"
rows <- c("beta", "uniform", "equidistant")

p11 <- tweak_fe(pick(boxplots, col1, rows[1]), title = col1_title, ic_mechanism = "beta", show_x = FALSE)
p12 <- tweak_fe(pick(boxplots, col2, rows[1]), title = col2_title, ic_mechanism = "beta", show_x = FALSE, show_y = FALSE)
p21 <- tweak_fe(pick(boxplots, col1, rows[2]), title = NULL, ic_mechanism = "uniform", show_x = FALSE)
p22 <- tweak_fe(pick(boxplots, col2, rows[2]), title = NULL, ic_mechanism = "uniform", show_x = FALSE, show_y = FALSE)
p31 <- tweak_fe(pick(boxplots, col1, rows[3]), title = NULL, ic_mechanism = "equidistant", show_x = TRUE)
p32 <- tweak_fe(pick(boxplots, col2, rows[3]), title = NULL, ic_mechanism = "equidistant", show_x = TRUE, show_y = FALSE)

ic_fe_boxplots <-
  (p11 | p12) /
  (p21 | p22) /
  (p31 | p32) +
  plot_layout(
    guides = "collect",
    design = "
      AA
      BB
      CC
    "
  ) &
  theme(
    legend.position      = "bottom",
    legend.direction     = "horizontal",
    legend.justification = "center",
    plot.margin          = margin(5, 5, 20, 5)
  )

png(file.path(dir_figures_fe, paste0("ic_fe_boxplots.png")), width = 12, height = 12, units = "in", res = 300)
print(ic_fe_boxplots)
dev.off()

### coverage barplots ----
barplots <- plot_coverage_fe(data = summary_fe, grouping_vars = c("problem", "ic_mechanism"))

if(!dir.exists(file.path(dir_figures, "coverages"))) {
  dir.create(file.path(dir_figures, "coverages"), recursive = TRUE)
}
for(i in 1:length(barplots)) {
  plot_name <- names(barplots)[i]

  ggsave(filename = paste0(dir_figures_fe, "/coverages/", plot_name, ".png"),
         plot = barplots[[i]], width = 10, height = 8)
}

## tables ----

### bias fe table ----
bias_summary_df <- create_bias_summary_fe(summary_fe)
bias_latex <- generate_bias_latex_fe(bias_summary_df)
writeLines(bias_latex, file.path(dir_tables, "ic-fe-bias.tex"))

### coverage fe table ----
coverage_summary_df <- create_coverage_summary_fe(summary_fe)
coverage_latex <- generate_coverage_latex_fe(coverage_summary_df)
writeLines(coverage_latex, file.path(dir_tables, "ic-fe-coverage.tex"))
