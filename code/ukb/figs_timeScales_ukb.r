
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

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb/baseline_plots"
dir_figures <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb/figures/export"

# Set font size for all plot elements
font_size <- 20
arrow <- " \u2192 "

# 0->1 and 0->4 plot ----
plots_pam_ssts_01 <- readRDS(file.path(dir_data, paste0("pam_ssts_plots_01.rds")))
plots_pam_ssts_04 <- readRDS(file.path(dir_data, paste0("pam_ssts_plots_04.rds")))
plots_pam_mts_01 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_01.rds")))
plots_pam_mts_04 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_04.rds")))

# --- 1. Prepare data for plotting ---

# Extract data from each ggplot object and combine
# Add columns to identify the model and transition, and format the transition label with an arrow
df_loghazard <- bind_rows(
  plots_pam_ssts_01$p_loghazard$data %>% mutate(model = "SSTS PAM", transition = "Healthy -> Mild CKD"),
  plots_pam_ssts_04$p_loghazard$data %>% mutate(model = "SSTS PAM", transition = "Healthy -> Death"),
  plots_pam_mts_01$p_loghazard$data %>% mutate(model = "MTS PAM", transition = "Healthy -> Mild CKD"),
  plots_pam_mts_04$p_loghazard$data %>% mutate(model = "MTS PAM", transition = "Healthy -> Death")
) %>%
  mutate(
    transition = gsub(" -> ", arrow, transition),
    model = factor(model, levels = c("SSTS PAM", "MTS PAM")))

df_tp <- bind_rows(
  plots_pam_ssts_01$p_tp$data %>% mutate(model = "SSTS PAM", transition = "Healthy -> Mild CKD"),
  plots_pam_ssts_04$p_tp$data %>% mutate(model = "SSTS PAM", transition = "Healthy -> Death"),
  plots_pam_mts_01$p_tp$data %>% mutate(model = "MTS PAM", transition = "Healthy -> Mild CKD"),
  plots_pam_mts_04$p_tp$data %>% mutate(model = "MTS PAM", transition = "Healthy -> Death")
) %>%
  mutate(
    transition = gsub(" -> ", arrow, transition),
    model = factor(model, levels = c("SSTS PAM", "MTS PAM")))

# --- 2. Create the two panel plots ---

# Define the new labels for the legend
transition_labels <- c("Healthy \u2192 Mild CKD", "Healthy \u2192 Death")
color_values <- c("blue", "darkgrey")
names(color_values) <- transition_labels

linetype_values <- c("longdash", "dotted")
names(linetype_values) <- c("SSTS PAM", "MTS PAM")

# Left Panel: Log-hazard
p_left <- ggplot(df_loghazard, aes(x = tend, y = loghazard, color = transition, linetype = model, fill = transition)) +
  geom_ribbon(aes(ymin = loghazard_lower, ymax = loghazard_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_values, breaks = transition_labels) +
  scale_fill_manual(values = color_values, breaks = transition_labels) +
  scale_linetype_manual(values = linetype_values) +
  labs(
    x = "Age",
    y = "Log-hazard",
    color = "Transition",
    linetype = "Model"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size),
    legend.key.width = unit(1.5, "cm")
  )

# Right Panel: Transition Probability
p_right <- ggplot(df_tp, aes(x = tend, y = trans_prob, color = transition, linetype = model, fill = transition)) +
  geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2, colour = NA, show.legend = FALSE) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = color_values, breaks = transition_labels) +
  scale_fill_manual(values = color_values, breaks = transition_labels) +
  scale_linetype_manual(values = linetype_values) +
  labs(
    x = "Age",
    y = "Transition Probability",
    color = "Transition",
    linetype = "Model"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = font_size),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size),
    legend.key.width = unit(1.5, "cm")
  )

# --- 3. Combine plots with patchwork ---

combined_plot <- p_left / p_right +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -50)
  )

ggsave(
  file.path(dir_figures, "ukb_baseline_01_04.png"),
  plot = combined_plot,
  width = 7,
  height = 12,
  units = "in",
  dpi = 300
)

# 1->2 and 1->4 plots ----

plots_pam_ssts_age_12 <- readRDS(file.path(dir_data, "pam_ssts_plots_age_12.rds"))
plots_pam_ssts_age_14 <- readRDS(file.path(dir_data, "pam_ssts_plots_age_14.rds"))
plots_pam_mts_age_12  <- readRDS(file.path(dir_data, "pam_mts_plots_age_12.rds"))
plots_pam_mts_age_14  <- readRDS(file.path(dir_data, "pam_mts_plots_age_14.rds"))

## contour plots
p12_ssts_log <- plots_pam_ssts_age_12$p_contour_loghazard
p12_ssts_tp  <- plots_pam_ssts_age_12$p_contour_tp
p12_mts_log  <- plots_pam_mts_age_12$p_contour_loghazard
p12_mts_tp   <- plots_pam_mts_age_12$p_contour_tp

p14_ssts_log <- plots_pam_ssts_age_14$p_contour_loghazard
p14_ssts_tp  <- plots_pam_ssts_age_14$p_contour_tp
p14_mts_log  <- plots_pam_mts_age_14$p_contour_loghazard
p14_mts_tp   <- plots_pam_mts_age_14$p_contour_tp

plot_list <- list(
  p12_tl = p12_ssts_log, p12_tr = p12_mts_log,
  p12_bl = p12_ssts_tp,  p12_br = p12_mts_tp,
  p14_tl = p14_ssts_log, p14_tr = p14_mts_log,
  p14_bl = p14_ssts_tp,  p14_br = p14_mts_tp
)

plot_list <- lapply(plot_list, function(p) {
  p + theme(
    text = element_text(size = font_size),
    plot.title = element_text(size = font_size, hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_blank(),
    legend.text = element_text(size = font_size),
    aspect.ratio = 1
  )
})

# --- Define plots with vertically stretched legends on the right ---
plot_list$p12_tl <- plot_list$p12_tl + labs(title = "SSTS PAM", subtitle = "Log-hazard", y = "Age at CKD onset", x = "Age") + theme(legend.position = "none")
plot_list$p12_tr <- plot_list$p12_tr + labs(title = "MTS PAM", subtitle = "Log-hazard", y = NULL, x = "Age") + theme(axis.text.y = element_blank(), legend.key.height = unit(1.8, "cm"))

plot_list$p12_bl <- plot_list$p12_bl + labs(title = NULL, subtitle = "Transition Probability", y = "Age at CKD onset", x = "Age") + theme(legend.position = "none")
plot_list$p12_br <- plot_list$p12_br + labs(title = NULL, subtitle = "Transition Probability", y = NULL, x = "Age") + theme(axis.text.y = element_blank(), legend.key.height = unit(1.8, "cm"))

plot_list$p14_tl <- plot_list$p14_tl + labs(title = "SSTS PAM", subtitle = "Log-hazard", y = "Age at CKD onset", x = "Age") + theme(legend.position = "none")
plot_list$p14_tr <- plot_list$p14_tr + labs(title = "MTS PAM", subtitle = "Log-hazard", y = NULL, x = "Age") + theme(axis.text.y = element_blank(), legend.key.height = unit(1.8, "cm"))

plot_list$p14_bl <- plot_list$p14_bl + labs(title = NULL, subtitle = "Transition Probability", y = "Age at CKD onset", x = "Age") + theme(legend.position = "none")
plot_list$p14_br <- plot_list$p14_br + labs(title = NULL, subtitle = "Transition Probability", y = NULL, x = "Age") + theme(axis.text.y = element_blank(), legend.key.height = unit(1.8, "cm"))

# --- Assemble the final panels from the completed rows ---
upper_panel <- (plot_list$p12_tl + plot_list$p12_tr) /
               (plot_list$p12_bl + plot_list$p12_br)

lower_panel <- (plot_list$p14_tl + plot_list$p14_tr) /
               (plot_list$p14_bl + plot_list$p14_br)

# --- Combine the final panels with their main titles ---
final_plot_12_14 <- wrap_elements(panel = upper_panel +
                                    plot_annotation(title = paste0("Transition Mild CKD", arrow, "Severe CKD"),
                                                    theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold")))) /
                    wrap_elements(panel = lower_panel +
                                    plot_annotation(title = paste0("Transition Mild CKD", arrow, "Death"),
                                                    theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold"))))

ggsave(
  file.path(dir_figures, "ukb_baseline_12_14_contour.png"),
  plot = final_plot_12_14,
  width = 10,
  height = 20,
  units = "in",
  dpi = 300
)

## slice plots
p12_ssts_log_slice <- plots_pam_ssts_age_12$p_slice_loghazard
p12_ssts_tp_slice  <- plots_pam_ssts_age_12$p_slice_tp
p12_mts_log_slice  <- plots_pam_mts_age_12$p_slice_loghazard
p12_mts_tp_slice   <- plots_pam_mts_age_12$p_slice_tp

p14_ssts_log_slice <- plots_pam_ssts_age_14$p_slice_loghazard
p14_ssts_tp_slice  <- plots_pam_ssts_age_14$p_slice_tp
p14_mts_log_slice  <- plots_pam_mts_age_14$p_slice_loghazard
p14_mts_tp_slice   <- plots_pam_mts_age_14$p_slice_tp

slice_plot_list <- list(
  p12_tl = p12_ssts_log_slice, p12_tr = p12_mts_log_slice,
  p12_bl = p12_ssts_tp_slice,  p12_br = p12_mts_tp_slice,
  p14_tl = p14_ssts_log_slice, p14_tr = p14_mts_log_slice,
  p14_bl = p14_ssts_tp_slice,  p14_br = p14_mts_tp_slice
)

slice_plot_list <- lapply(slice_plot_list, function(p) {
  p + theme(
    text = element_text(size = font_size),
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_text(size = font_size),
    axis.text = element_text(size = font_size),
    legend.title = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    legend.position = "none"
  )
})

slice_plot_list$p12_tl <- slice_plot_list$p12_tl + labs(title = "SSTS PAM", x = NULL) + theme(plot.title = element_text(size = font_size, hjust = 0.5))
slice_plot_list$p12_tr <- slice_plot_list$p12_tr + labs(title = "MTS PAM", y = NULL, x = NULL) + theme(axis.text.y = element_blank(), plot.title = element_text(size = font_size, hjust = 0.5))
slice_plot_list$p12_bl <- slice_plot_list$p12_bl + labs(title = NULL, x = "Age", y = "Transition\nProbability")
slice_plot_list$p12_br <- slice_plot_list$p12_br + labs(title = NULL, y = NULL, x = "Age") + theme(axis.text.y = element_blank())

slice_plot_list$p14_tl <- slice_plot_list$p14_tl + labs(title = "SSTS PAM", x = NULL) + theme(plot.title = element_text(size = font_size, hjust = 0.5))
slice_plot_list$p14_tr <- slice_plot_list$p14_tr + labs(title = "MTS PAM", y = NULL, x = NULL) + theme(axis.text.y = element_blank(), plot.title = element_text(size = font_size, hjust = 0.5))
slice_plot_list$p14_bl <- slice_plot_list$p14_bl + labs(title = NULL, x = "Age", y = "Transition\nProbability")
slice_plot_list$p14_br <- slice_plot_list$p14_br + labs(title = NULL, y = NULL, x = "Age") + theme(axis.text.y = element_blank())

upper_panel_slice <- (slice_plot_list$p12_tl + slice_plot_list$p12_tr) /
                     (slice_plot_list$p12_bl + slice_plot_list$p12_br) +
                     plot_layout(guides = "collect") & theme(legend.position = "bottom")

lower_panel_slice <- (slice_plot_list$p14_tl + slice_plot_list$p14_tr) /
                     (slice_plot_list$p14_bl + slice_plot_list$p14_br) +
                     plot_layout(guides = "collect") & theme(legend.position = "bottom")

final_plot_12_14_slice <- wrap_elements(panel = upper_panel_slice +
                                          plot_annotation(title = paste0("Transition Mild CKD", arrow, "Severe CKD"),
                                                          theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold")))) /
                          wrap_elements(panel = lower_panel_slice +
                                          plot_annotation(title = paste0("Transition Mild CKD", arrow, "Death"),
                                                          theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold"))))

ggsave(
  file.path(dir_figures, "ukb_baseline_12_14_slice.png"),
  plot = final_plot_12_14_slice,
  width = 10,
  height = 12,
  units = "in",
  dpi = 300
)

# 2->3 and 2->4 plot ----

## contour plots

for (trans_num in c("23", "24")) {

  ssts_file <- file.path(dir_data, paste0("pam_ssts_plots_age_", trans_num, ".rds"))
  mts_file  <- file.path(dir_data, paste0("pam_mts_plots_age_", trans_num, ".rds"))

  plots_pam_ssts <- readRDS(ssts_file)
  plots_pam_mts  <- readRDS(mts_file)

  p_ssts_log <- plots_pam_ssts$p_contour_loghazard
  p_ssts_tp  <- plots_pam_ssts$p_contour_tp
  p_mts_log  <- plots_pam_mts$p_contour_loghazard
  p_mts_tp   <- plots_pam_mts$p_contour_tp

  p_tl <- p_ssts_log +
    labs(title = "SSTS PAM", subtitle = "Log-hazard", y = "Age at CKD progression", x = "Age") +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.subtitle = element_text(size = font_size, face = "italic"),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      legend.position = "none",
      aspect.ratio = 1,
      plot.margin = unit(c(5.5, 2, 5.5, 1), "pt")
    )

  p_tr <- p_mts_log +
    labs(title = "MTS PAM", subtitle = "Log-hazard", y = NULL, x = "Age") +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.subtitle = element_text(size = font_size, face = "italic"),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1,
      plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")
    )

  p_bl <- p_ssts_tp +
    labs(title = "SSTS PAM", subtitle = "Transition Probability", y = "Age at CKD progression", x = "Age") +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.subtitle = element_text(size = font_size, face = "italic"),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      legend.position = "none",
      aspect.ratio = 1,
      plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")
    )

  p_br <- p_mts_tp +
    labs(title = "MTS PAM", subtitle = "Transition Probability", y = NULL, x = "Age") +
    theme(
      text = element_text(size = font_size),
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.subtitle = element_text(size = font_size, face = "italic"),
      axis.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      legend.text = element_text(size = font_size),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      aspect.ratio = 1,
      plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt")
    )

  top_row_unit <- (p_tl + p_tr) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.key.width = unit(2.5, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = font_size)
    )

  bottom_row_unit <- (p_bl + p_br) +
    plot_layout(guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.key.width = unit(2.5, "cm"),
      legend.title = element_blank(),
      legend.text = element_text(size = font_size)
    )

  combined_plot <- top_row_unit / bottom_row_unit

  from_digit <- substr(trans_num, 1, 1)
  from_state <- ifelse(from_digit == "2", "Severe CKD", "Other")
  to_digit   <- substr(trans_num, 2, 2)
  to_state   <- ifelse(to_digit == "3", "ESKD", "Death")
  main_title <- paste0("Transition ", from_state, arrow, to_state)

  final_plot <- combined_plot +
    plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold"))
    )

  output_filename <- paste0("ukb_baseline_", trans_num, "_contour.png")
  ggsave(
    file.path(dir_figures, output_filename),
    plot = final_plot,
    width = 16,
    height = 20,
    units = "in",
    dpi = 300
  )
}

## slice plots

for (trans_num in c("23", "24")) {

  ssts_file <- file.path(dir_data, paste0("pam_ssts_plots_age_", trans_num, ".rds"))
  mts_file  <- file.path(dir_data, paste0("pam_mts_plots_age_", trans_num, ".rds"))

  plots_pam_ssts <- readRDS(ssts_file)
  plots_pam_mts  <- readRDS(mts_file)

  p_ssts_log <- plots_pam_ssts$p_slice_loghazard
  p_ssts_tp  <- plots_pam_ssts$p_slice_tp
  p_mts_log  <- plots_pam_mts$p_slice_loghazard
  p_mts_tp   <- plots_pam_mts$p_slice_tp

  plot_list <- list(
    p_tl = p_ssts_log, p_tr = p_mts_log,
    p_bl = p_ssts_tp,  p_br = p_mts_tp
  )

  plot_list <- lapply(plot_list, function(p) {
    p + theme(
        text = element_text(size = font_size),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.title = element_text(size = font_size),
        axis.text = element_text(size = font_size),
        legend.title = element_text(size = font_size),
        legend.text = element_text(size = font_size),
        legend.position = "none"
      )
  })

  plot_list$p_tl <- plot_list$p_tl +
    labs(title = "SSTS PAM", x = NULL) +
    theme(
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.margin = unit(c(5.5, 5.5, 20, 5.5), "pt")
    )

  plot_list$p_tr <- plot_list$p_tr +
    labs(title = "MTS PAM", y = NULL, x = NULL) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = font_size, hjust = 0.5),
      plot.margin = unit(c(5.5, 5.5, 20, 5.5), "pt")
    )

  plot_list$p_bl <- plot_list$p_bl +
    labs(title = "SSTS PAM", x = "Age", y = "Transition Probability") +
    theme(plot.title = element_text(size = font_size, hjust = 0.5))


  plot_list$p_br <- plot_list$p_br +
    labs(title = "MTS PAM", y = NULL, x = "Age") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(size = font_size, hjust = 0.5)
    )

  combined_plot <- (plot_list$p_tl + plot_list$p_tr) / (plot_list$p_bl + plot_list$p_br) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  from_digit <- substr(trans_num, 1, 1)
  from_state <- ifelse(from_digit == "2", "Severe CKD", "Other")
  to_digit   <- substr(trans_num, 2, 2)
  to_state   <- ifelse(to_digit == "3", "ESKD", "Death")
  main_title <- paste0("Transition ", from_state, arrow, to_state)

  final_plot <- combined_plot +
    plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = font_size + 2, hjust = 0.5, face = "bold"))
    )

  output_filename <- paste0("ukb_baseline_", trans_num, "_slice.png")
  ggsave(
    file.path(dir_figures, output_filename),
    plot = final_plot,
    width = 16,
    height = 14,
    units = "in",
    dpi = 300
  )
}
