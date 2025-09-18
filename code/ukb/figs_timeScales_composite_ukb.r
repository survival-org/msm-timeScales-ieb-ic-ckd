
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/ukb/helpers_ukb.r")

library(etm)
library(mvna)
library(kmi)
library(mgcv)
library(scales)
library(nnet)
library(patchwork)

# set directories ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_data_baseline_plots <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb/baseline_plots/composite"
dir_figures <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb/figures/export"

# Set font size for all plot elements
font_size <- 32
arrow <- " \u2192 "

# create individual plots ----

## transition probabilities ----
ped_new <- readRDS(file.path(dir_data, "ped_pam_ssts.rds"))

age_onset_slices  <- c(50, 60)
prog_after_onset_slices <- c(2, 5, 10, 15)

ranges <- extract_ranges(ped_new, ped_new, transitions = c("0->1", "1->2", "2->3"), scales = c("loghazard", "hazard", "cumu_hazard"))

plot_tp_01 <- create_2d_plots(ped_new, model, "0->1", ranges, save_plots = FALSE)$p_tp
plot_tp_12 <- create_3d_plots(ped_new, model, "1->2", ranges, "age", age_onset_slices, prog_after_onset_slices, save_plots = FALSE)$p_slice_tp
plot_tp_23 <- create_3d_plots(ped_new, model, "2->3", ranges, "age", age_onset_slices, prog_after_onset_slices, save_plots = FALSE)$p_slice_tp

p_tp_01 <- ggplot(plot_tp_01$data, aes(x = tend, y = trans_prob)) +
  geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2, color = NA, fill = "blue", show.legend = FALSE) +
  geom_line(color = "blue", linewidth = 1.2) +
  # scale_color_manual(values = "blue", breaks = "Healthy \u2192 Mild CKD") +
  # scale_fill_manual(values = "blue", breaks = "Healthy \u2192 Mild CKD") +
  scale_linetype_manual(values = "solid") +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "Healthy \u2192 Mild CKD",
    x = "Age",
    y = "Transition Probability"
  ) +
  theme_bw(base_size = font_size) +
  theme(
    text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size),
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = font_size, hjust = 0.5),
    plot.margin = unit(c(0,0,0,0), "pt")
  )

full_color_palette <- brewer.pal(length(c(40, 50, 60, 70)), "RdBu")
color_map <- c(
  "50" = full_color_palette[2], # The light orange
  "60" = full_color_palette[3]  # The light blue
)

p_tp_12 <- plot_tp_12 +
  scale_colour_manual(
    name = "Age at CKD onset", # Make sure the legend title is correct
    values = color_map
  ) +
  scale_fill_manual(
    name = "Age at CKD onset", # Match the legend title
    values = color_map
  ) +
  labs(
    title = "Mild CKD \u2192 Severe CKD",
    x = "Age",
    y = "Transition Probability") +
  theme_bw(base_size = font_size) +
  theme(
    text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.position = "bottom",
    plot.title = element_text(size = font_size, hjust = 0.5),
    plot.margin = unit(c(0,0,0,0), "pt")
  )

p_tp_23 <- plot_tp_23 +
  labs(
    title = "Severe CKD \u2192 ESKD",
    x = "Age",
    y = "Transition Probability") +
  theme_bw(base_size = font_size) +
  theme(
    text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.position = "bottom",
    plot.title = element_text(size = font_size, hjust = 0.5),
    plot.margin = unit(c(0,0,0,0), "pt")
  )

## left truncation ----
ssts_onset_df <- readRDS(file.path(dir_data, "leftTruncation_ssts_onset_df.rds")) %>%
  filter(transition_after_onset_strat == "1->2")
ssts_prog_df <- readRDS(file.path(dir_data, "leftTruncation_ssts_prog_df.rds")) %>%
  filter(transition_after_progression == "eskd")

plot_leftTrunc_12 <- ggplot(
  ssts_onset_df, aes(x = .data[["age_onset"]], y = .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "#0072B2", colour = NA) +
    geom_line(color = "#0072B2", linewidth = 1.2) +
    scale_x_continuous(
      limits = c(35, NA),
      breaks = c(40, 50, 60, 70, 80),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    scale_y_continuous(
      limits = c(-4, 4),
      breaks = seq(-4, 4, by = 2)
    ) +
    labs(
      title = "Mild CKD \u2192 Severe CKD",
      x = "Age at CKD Onset",
      y = "Log-hazard"
    ) +
    theme_bw(base_size = font_size) +
    theme(
      text = element_text(size = font_size),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      strip.text = element_text(size = font_size, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      legend.position = "none",
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.margin = unit(c(0,0,0,0), "pt")
    )

plot_leftTrunc_23 <- ggplot(
  ssts_prog_df, aes(x = .data[["age_progression"]], y = .estimate)) +
    geom_ribbon(aes(ymin = .lower_ci, ymax = .upper_ci), alpha = 0.2, fill = "#0072B2", colour = NA) +
    geom_line(color = "#0072B2", linewidth = 1.2) +
    scale_x_continuous(
      limits = c(35, NA),
      breaks = c(40, 50, 60, 70, 80),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    scale_y_continuous(
      limits = c(-4, 4),
      breaks = seq(-4, 4, by = 2)
    ) +
    labs(
      title = "Severe CKD \u2192 ESKD",
      x = "Age at CKD Progression",
      y = "Log-hazard"
    ) +
    theme_bw(base_size = font_size) +
    theme(
      text = element_text(size = font_size),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      strip.text = element_text(size = font_size, face = "bold"),
      strip.background = element_rect(fill = "grey90", color = NA),
      legend.position = "none",
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.margin = unit(c(0,0,0,0), "pt")
    )

# final panel plot ----

## transition probability panel
panel_transProb <- (p_tp_01 + p_tp_12) / p_tp_23 +
  plot_annotation(
    title = 'Transition Probabilities',
    theme = theme(plot.title = element_text(size = font_size * 1.2, hjust = 0.5, face = "bold"))
  )

# ggsave(
#   file.path(dir_figures, "ukb_composite_transProb.png"),
#   plot = panel_transProb,
#   width = 16, height = 14, units = "in", dpi = 300
# )

## left truncation panel
panel_leftTrunc <- (plot_leftTrunc_12 + plot_leftTrunc_23) +
  plot_annotation(
    title = 'Centered Effects of Left-truncation Times',
    theme = theme(plot.title = element_text(size = font_size * 1.2, hjust = 0.5, face = "bold"))
  )

# ggsave(
#   file.path(dir_figures, "ukb_composite_leftTrunc.png"),
#   plot = panel_leftTrunc,
#   width = 16, height = 10, units = "in", dpi = 300
# )

## combined
panel_combined <- wrap_plots(
  panel_transProb,
  panel_leftTrunc,
  ncol = 1
) +
  plot_layout(heights = c(2.5, 1)) +
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(face = "bold", size = font_size + 4))

ggsave(
  file.path(dir_figures, "ukb_composite_transProb_leftTrunc.png"),
  plot = panel_combined,
  width = 22, height = 24, units = "in", dpi = 300
)
