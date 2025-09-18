
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/ukb/helpers_ukb.r")

library(mgcv)
library(ggplot2)
library(gratia)
library(stringr)


# set directories ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_figures <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb/figures/export"

# load data ----
# min_age <- 35

# ped_events <- readRDS(file.path(dir_data, file_ped_events)) %>%
#   filter(tstart >= min_age) %>%
#   mutate(
#     to = ifelse(to == "death", 4, to),
#     transition = paste0(from, "->", to) %>% as.factor(), # for pammtools::add_surv_prob
#     diabetes = ifelse(diabetes == 0, 0, 1)
#   ) %>%
#   add_transVars()

pam_ssts <- readRDS(file.path(dir_data, "pam_ssts_ukb.rds"))
pam_mts <- readRDS(file.path(dir_data, "pam_mts_ukb.rds"))

# fit models ----
# formula_ssts <- "ped_status ~
#   s(tend, bs = 'ps', k = 20, by = transition) +
#   s(age_onset, bs = 'ps', k = 20, by = transition_after_onset_strat) +
#   s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
#   sex*transition"

# pam_ssts <- mgcv::bam(
#   formula = as.formula(formula_ssts),
#   data = ped_events,
#   family = poisson(),
#   offset = offset,
#   method = "fREML",
#   discrete = TRUE)

# saveRDS(pam_ssts, file = file.path(dir_data, "pam_ssts_ukb.rds"))

# ### pam mts ----
# formula_mts <- "ped_status ~
#   s(tend, bs = 'ps', k = 20, by = transition_to_death) +
#   s(tend_onset, bs = 'ps', k = 20, by = transition_after_onset) +
#   s(age_onset, bs = 'ps', k = 20, by = transition_after_onset) +
#   s(tend_progression, bs = 'ps', k = 20, by = transition_after_progression) +
#   s(age_progression, bs = 'ps', k = 20, by = transition_after_progression) +
#   sex*transition"

# pam_mts <- mgcv::bam(
#   formula = as.formula(formula_mts),
#   data = ped_events,
#   family = poisson(),
#   offset = offset,
#   method = "fREML",
#   discrete = TRUE)

# saveRDS(pam_mts, file = file.path(dir_data, "pam_mts_ukb.rds"))

# create individual plots ----

font_size <- 20
y_axis_label <- "Centered Effect (Log-hazard Scale)"

new_color_map <- c(
  "Mild CKD \u2192 Severe CKD" = "#0072B2",
  "Mild CKD \u2192 Death"      = "#5A5A5A",
  "Severe CKD \u2192 ESKD"     = "#0072B2",
  "Severe CKD \u2192 Death"     = "#5A5A5A"
)

extract_and_process <- function(model, covariate_name, transition_mapping) {
  smooth_labels <- sapply(model$smooth, function(s) s$label)

  all_data <- lapply(names(transition_mapping), function(final_transition_name) {
    search_pattern <- transition_mapping[[final_transition_name]]

    target_smooth <- smooth_labels[
      str_detect(smooth_labels, covariate_name) & str_detect(smooth_labels, search_pattern)
    ]

    if (length(target_smooth) == 0) return(NULL)

    smooth_estimates(model, select = target_smooth) %>%
      mutate(
        .lower_ci = .estimate - (1.96 * .se),
        .upper_ci = .estimate + (1.96 * .se),
        transition = final_transition_name
      )
  })

  bind_rows(all_data)
}

ssts_onset_map <- c("Mild CKD \u2192 Severe CKD" = "1->2", "Mild CKD \u2192 Death" = "1->4")
mts_onset_map  <- c("Mild CKD \u2192 Severe CKD" = "progression", "Mild CKD \u2192 Death" = "death")
ssts_prog_map  <- c("Severe CKD \u2192 ESKD" = "eskd", "Severe CKD \u2192 Death" = "death")
mts_prog_map   <- c("Severe CKD \u2192 ESKD" = "eskd", "Severe CKD \u2192 Death" = "death")

ssts_onset_df <- extract_and_process(pam_ssts, "age_onset", ssts_onset_map)
mts_onset_df  <- extract_and_process(pam_mts, "age_onset", mts_onset_map)
ssts_prog_df  <- extract_and_process(pam_ssts, "age_progression", ssts_prog_map)
mts_prog_df   <- extract_and_process(pam_mts, "age_progression", mts_prog_map)

saveRDS(ssts_onset_df, file.path(dir_data, "leftTruncation_ssts_onset_df.rds"))
saveRDS(mts_onset_df, file.path(dir_data, "leftTruncation_mts_onset_df.rds"))
saveRDS(ssts_prog_df, file.path(dir_data, "leftTruncation_ssts_prog_df.rds"))
saveRDS(mts_prog_df, file.path(dir_data, "leftTruncation_mts_prog_df.rds"))

create_facet_plot <- function(data, covariate_col, plot_title, x_axis_label, transition_order) {
  if (nrow(data) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(plot_title, "\n(No data found)")))
  }

  data$transition <- factor(data$transition, levels = transition_order)

  ggplot(data, aes(x = .data[[covariate_col]], y = .estimate, color = transition, fill = transition)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1.2) +
    facet_wrap(~ transition, scales = "free_y", ncol = 1) +
    scale_color_manual(values = new_color_map) +
    scale_fill_manual(values = new_color_map) +
    scale_x_continuous(
      limits = c(35, NA),
      breaks = c(40, 50, 60, 70, 80),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    labs(
      title = plot_title,
      x = x_axis_label,
      y = y_axis_label
    ) +
    theme_bw(base_size = font_size) +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(hjust = 0.5, face = "bold", size = font_size + 2),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      strip.text = element_text(size = font_size, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      legend.position = "none"
    )
}

plot_ssts_onset <- create_facet_plot(
  data = ssts_onset_df,
  covariate_col = "age_onset",
  plot_title = "SSTS PAM",
  x_axis_label = "Age at CKD Onset",
  transition_order = c("Mild CKD \u2192 Severe CKD", "Mild CKD \u2192 Death")
)
ggsave(
  file.path(dir_figures, "ukb_leftTrunc_ssts_onset.png"),
  plot = plot_ssts_onset,
  width = 8, height = 10, units = "in", dpi = 300
)

plot_mts_onset <- create_facet_plot(
  data = mts_onset_df,
  covariate_col = "age_onset",
  plot_title = "MTS PAM",
  x_axis_label = "Age at CKD Onset",
  transition_order = c("Mild CKD \u2192 Severe CKD", "Mild CKD \u2192 Death")
)
ggsave(
  file.path(dir_figures, "ukb_leftTrunc_mts_onset.png"),
  plot = plot_mts_onset,
  width = 8, height = 10, units = "in", dpi = 300
)

plot_ssts_prog <- create_facet_plot(
  data = ssts_prog_df,
  covariate_col = "age_progression",
  plot_title = "SSTS PAM",
  x_axis_label = "Age at CKD Progression",
  transition_order = c("Severe CKD \u2192 ESKD", "Severe CKD \u2192 Death")
)
ggsave(
  file.path(dir_figures, "ukb_leftTrunc_ssts_prog.png"),
  plot = plot_ssts_prog,
  width = 8, height = 10, units = "in", dpi = 300
)

plot_mts_prog <- create_facet_plot(
  data = mts_prog_df,
  covariate_col = "age_progression",
  plot_title = "MTS PAM",
  x_axis_label = "Age at CKD Progression",
  transition_order = c("Severe CKD \u2192 ESKD", "Severe CKD \u2192 Death")
)
ggsave(
  file.path(dir_figures, "ukb_leftTrunc_mts_prog.png"),
  plot = plot_mts_prog,
  width = 8, height = 10, units = "in", dpi = 300
)