library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)
library(patchwork)

# setup ----
setwd("nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd")
source("code/simulations/helpers_sim.r")
source("code/simulations/helpers_ieb.r")

registry <- "results/simulations/registries/sim-ieb-registry"
dir_datasets <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/datasets/"
dir_figures <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/figures/ieb/"
dir_tables <- "/nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd/results/simulations/tables/ieb/"

# load ----
res_msm <- readRDS(paste0(dir_datasets, "sim-ieb-results_msm.rds"))
res_cor <- readRDS(paste0(dir_datasets, "sim-ieb-results_cor.rds"))

# plots ----
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
      plot   = p1[[1]] / p1[[2]],
      width  = 8, height = 8
    )

  #   # 2) coverage barplot
  #   p2 <- plot_coverage(df_sub)
  #   ggsave(
  #     filename = file.path(
  #       dir_figures,
  #       "coverages",
  #       paste0("coverage_boxplots_",
  #              dist, "_", scn, ".png")
  #     ),
  #     plot   = p2,
  #     width  = 8, height = 8
  #   )
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
  labs(title = expression(atop("Moderate Effect Sizes",
                               x[1] ~ " ~ Ber(0.5) & " ~ x[2] ~ " ~ Ber(0.5)")))

# 3. Filter data for the RIGHT plot
df_right <- res_cor %>%
  filter(dist_x1_x2    == "normal1_normal1",
         beta_scenario == "extreme")

# 4. Create the RIGHT plot and add the plotmath title
p_right <- plot_cor(df_right, font_size = 20) +
  labs(title = expression(atop("Very Strong Effect Sizes",
                               x[1] ~ " ~ N(0,1) & " ~ x[2] ~ " ~ N(0,1)"))) +
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
  filename = file.path(dir_figures, "ieb_cor.png"),
  plot   = combined_plot, width = 12, height = 6, dpi = 300
)

### coefficient facet boxplot ----

# 1. Get the LISTS of plots from your function
# This returns a list of 2 ggplot objects for each call
df_left <- res_msm %>%
  filter(dist_x1_x2    == "bernoulli0.5_bernoulli0.5",
         beta_scenario == "moderate")

df_right <- res_msm %>%
  filter(dist_x1_x2    == "normal1_normal1",
         beta_scenario == "strong")

p_list_left <- plot_coefs(df_left, font_size = 20)
p_list_right <- plot_coefs(df_right, font_size = 20)

ggsave("test.png", p_list_left[[1]], width = 8, height = 8, dpi = 300)

# 2. Add overall titles as separate "plots"
title_left_text <- expression(atop("Moderate Effect Sizes", x[1] ~ " ~ Ber(0.5) & " ~ x[2] ~ " ~ Ber(0.5)"))
title_left <- ggplot() + labs(title = title_left_text) + theme_void() +
              theme(plot.title = element_text(hjust = 0.5, size = 22))

title_right_text <- expression(atop("Very Strong Effect Sizes", x[1] ~ " ~ N(0,1) & " ~ x[2] ~ " ~ N(0,1)"))
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

## effect sizes ----

# UPDATE HERE!!! COPY FROM simulations_ieb.r
betas <- list(
  moderate = c(beta_0_01 = -3.9, beta_0_03 = -4.0,
               beta_0_12 = -3.4, beta_0_13 = -3.4,
               beta_1_01 = 0.2, beta_1_03 = 0.1,
               beta_1_12 = 0.2, beta_1_13 = 0.1,
               beta_2_01 = 0.2, beta_2_03 = 0.0,
               beta_2_12 = 0.2, beta_2_13 = 0.0),
  strong   = c(beta_0_01 = -3.9, beta_0_03 = -4.0,
               beta_0_12 = -3.4, beta_0_13 = -3.4,
               beta_1_01 = 0.4, beta_1_03 = 0.0,
               beta_1_12 = 0.4, beta_1_13 = 0.0,
               beta_2_01 = 0.4, beta_2_03 = 0.0,
               beta_2_12 = 0.4, beta_2_13 = 0.0),
  extreme  = c(beta_0_01 = -3.9, beta_0_03 = -4.0,
               beta_0_12 = -4.4, beta_0_13 = -3.4,
               beta_1_01 = 0.6, beta_1_03 = 0.0,
               beta_1_12 = 0.6, beta_1_13 = 0.0,
               beta_2_01 = 0.6, beta_2_03 = 0.0,
               beta_2_12 = 0.6, beta_2_13 = 0.0)
)

latex_effect_sizes <- convert_to_latex_effect_sizes(betas)
writeLines(latex_effect_sizes, file.path(dir_tables, "ieb-effect-sizes.tex"))

## correlation ----
cor_summary <- res_cor %>%
  group_by(problem, dist_x1_x2, beta_scenario, from) %>%
  summarise(
    rho_mean = mean(rho),
    rho_ci = list(round(quantile(rho, probs = c(0.025, 0.975)), 4))
    ) %>%
    mutate(
      rho_lower = map_dbl(rho_ci, ~ .[1]),
      rho_upper = map_dbl(rho_ci, ~ .[2]),
      rho_ci = map_chr(rho_ci, ~ paste0("(", .[1], ", ", .[2], ")")))

processed_df_cor <- process_correlation_summary(cor_summary)
latex_cor <- convert_to_latex_cor(processed_df_cor, colsep = "2pt")
writeLines(latex_cor, file.path(dir_tables, "ieb_cor.tex"))

## bias ----
ieb_coef_bias_cov_summary <- create_csv_summary(res_msm)
View(ieb_coef_bias_cov_summary)

# CHECK!!!
write.csv(ieb_coef_bias_cov_summary, file.path(dir_tables, "ieb_coef_bias_cov_summary.csv"), row.names = FALSE)
