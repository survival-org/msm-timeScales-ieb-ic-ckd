
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/helpers_ukb.r")

library(etm)
library(mvna)
library(kmi)
library(mgcv)
library(scales)
library(nnet)

# set directories ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_out <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb/msm"
file_pheno <- "pheno_ukb_noAKI6months_noNephrectomy6months_ni2+_GPclinical.txt"
file_events <- "events_ukb_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"
file_ped_events <- "ped_events_ukb_noAKI6months_noNephrectomy6months_ni2+_GPclinical_mid_imputed_r0_ageScale.rds"

# load data ----
pheno <- fread(file.path(dir_data, file_pheno))
min_age <- 35

events <- readRDS(file.path(dir_data, file_events)) %>%
  mutate(to = ifelse(to == "death", 4, to), transition = paste0(from, "->", to)) %>%
  filter(agestart >= min_age)

ped_events <- readRDS(file.path(dir_data, file_ped_events)) %>%
  filter(tstart >= min_age) %>%
  mutate(
    to = ifelse(to == "death", 4, to),
    transition = paste0(from, "->", to) %>% as.factor(), # for pammtools::add_surv_prob
  ) %>%
  add_transVars()

# cols_to_include <- c("id", "tstart", "tend", "interval", "offset", "ped_status", "from", "to", "transition", "sex",
#                      "tend_onset", "age_onset", "transition_after_onset", "transition_after_onset_strat",
#                      "tend_progression", "age_progression", "transition_after_progression")

# ped_events_anon <- ped_events %>%
#   select(all_of(cols_to_include))

# Descriptive statistics ----

## Table 1
dim(events)
table(events$transition)
idx <- unique(events$id)
length(idx)
events %>% filter(transition %in% c("0->1", "1->2", "2->3")) %>% nrow
events %>% filter(to == 4) %>% nrow
events %>% group_by(id) %>% slice(1) %>% ungroup() %>% pull(sex) %>% table()
events %>% group_by(id) %>% slice(1) %>% ungroup() %>% pull(agestart) %>% summary()
events %>% group_by(id) %>% slice(n()) %>% ungroup() %>% pull(agestop) %>% summary()

## Histogram of transition times ----

events <- events %>%
  mutate(
    t = (agestop - agestart - 1)/3.6
  )

for(trans in unique(events$transition)) {
  p <- ggplot(
    events %>% filter(transition == trans),
    aes(x = t)) +
    geom_histogram(
      bins = 30,
      fill = "steelblue",
      color = "black",
      boundary = 0
    ) +
    ggtitle(paste("Transition time histogram for", trans))

  ggsave(file.path(dir_out, paste0("/figures/descriptives/histogram_", trans, ".png")), plot = p, width = 10, height = 6, units = "cm")
}

## Stacked probability / proportion plot (DATAMAP guy from SAfJR, Friday morning), using mstate package ----

# 0) parameters & prep
global_min_age <- 30
global_max_age <- 80

# 0a) IDs with at least one non‐censoring event
ids_event <- unique(events[to != "cens", id])

# 1) load & subset events
dt_all <- as.data.table(events)[
  id %in% ids_event &               # only those with any event
  !(to == "cens" & from == "0")     # drop only 0->cens
]

# 1b) truncate truly early records
dt_all <- dt_all[agestart >= global_min_age]

# 2) build per‐id occupancy by tiling each interval
occ <- dt_all[, {
  tr         <- .SD[order(agestart)]
  rows       <- list()
  for (i in seq_len(nrow(tr))) {
    r <- tr[i]
    # tile [floor(agestart) … ceiling(agestop)-1] in the 'from' state
    a0 <- floor(r$agestart)
    a1 <- ceiling(r$agestop) - 1
    if (a0 <= a1) {
      rows[[length(rows)+1]] <- data.table(
        age   = a0:a1,
        state = as.character(r$from)
      )
    }
    # carry forward only if absorbing (3 or 4)
    if (r$to %in% c("3","4")) {
      carry0 <- ceiling(r$agestop)
      if (carry0 <= global_max_age) {
        rows[[length(rows)+1]] <- data.table(
          age   = carry0:global_max_age,
          state = as.character(r$to)
        )
      }
      break
    }
  }
  rbindlist(rows)
}, by = id]

# remove any duplicates
occ <- unique(occ, by = c("id","age"))

# 3) count how many ids in each (age, state)
counts <- occ[, .(n = uniqueN(id)), by = .(age, state)]

# 4) complete the grid & fill zeros
all_states <- as.character(0:4)
ages      <- seq(global_min_age, global_max_age)
grid      <- CJ(age = ages, state = all_states)

df_plot <- grid[counts, on=.(age, state)]
df_plot[is.na(n), n := 0L]

# 5) compute total, prob, prop
df_plot[, total := sum(n), by = age]
df_plot[, prob  := n / total]
df_plot[, prop  := 100 * prob]

# 6) factor state for consistent stacking
df_plot[, state := factor(state, levels = all_states)]

# 7) single labels at midpoint age
mid_age <- floor((global_min_age + global_max_age)/2)
label_df <- df_plot[age == mid_age]
label_df[, y_prob := cumsum(prob) - prob/2, by = state]
label_df[, y_prop := cumsum(prop) - prop/2, by = state]

# 8a) probability plot with x‐ticks at 30,40,…,80
p_prob <- ggplot(
  df_plot,
  aes(x = age, y = prob, fill = state)) +
  geom_area() +
  # geom_text(
  #   data    = label_df,
  #   aes(x = age, y = y_prob, label = state),
  #   color   = "white",
  #   size    = 4
  # ) +
  scale_x_continuous(limits = c(global_min_age, global_max_age), breaks = seq(30, 80, by=10)) +
  # scale_y_continuous(labels = percent_format(1)) +
  labs(
    x     = "Age",
    y     = "Probability",
    fill  = "State"
  )

p_prob <- ggplot(df_plot, aes(
    x    = age,
    y    = prob,
    fill = state
  )) +
  geom_col(
    position = "stack",
    width    = 1       # full–width bars for each integer age
  ) +
  scale_x_continuous(
    breaks = seq(30, 80, by = 10),
    limits = c(30, 80),
    expand = c(0,0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(1)
  ) +
  labs(x = "Age", y = "Probability", fill = "State")

ggsave(file.path(dir_out, "/figures/stacked_probability_plot.png"), plot = p_prob, width = 10, height = 6, units = "cm")
ggsave(file.path(dir_out, "/figures/stacked_proportions_plot.png"), plot = p_prop, width = 10, height = 6, units = "cm")

## distribution of risk factors ----

ped_baseline <- ped_events %>%
  group_by(id, from) %>%
  slice_min(order_by = date_at_exam, with_ties = FALSE) %>%  # pick the single earliest row
  ungroup() %>%
  select(id, from, sc_UACR, sc_BMI, smoking)

events_w <- events %>%
  # bring in the two new vars
  left_join(ped_baseline, by = c("id", "from")) %>%
  # make from a nice factor (optional labels)
  mutate(from = factor(from,
                       levels = c(0,1,2),
                       labels = c("State 0","State 1","State 2"))) %>%
  # compute 1/n per state for weighting
  group_by(from) %>%
  mutate(weight = 1 / n()) %>%
  ungroup() %>%
    mutate(
      albuminuria = case_when(
        sc_UACR >= 30 & sc_UACR < 300 ~ "micro",
        sc_UACR >= 300 ~ "macro",
        TRUE ~ "normo"
      ) %>%
      factor(levels = c("normo", "micro", "macro")),
    BMI_cat = case_when(
      sc_BMI >= 25 & sc_BMI < 30 ~ "overweight",
      sc_BMI >= 30 ~ "obese",
      TRUE ~ "normal") %>%
      factor(levels = c("normal", "overweight", "obese")),
    rs_77924615_A_G_genotype = case_when(
      rs77924615_A_G > 1.5 ~ 2,
      rs77924615_A_G > 0.5 ~ 1,
      TRUE ~ 0
    ) %>%
    factor(levels = c(0, 1, 2)),
    smoking = factor(smoking, levels = c(0, 1)),
    diabetes = factor(diabetes, levels = c(0, 1))
  )

events_w2 <- events_w %>%
  # compute weight_g = 1/(#obs in each state & genotype)
  group_by(from, rs_77924615_A_G_genotype) %>%
  mutate(weight_g = 1 / n()) %>%
  ungroup()

library(nnet)
library(splines)

events_ps <- events_w2 %>%
  group_by(from) %>%
  do({
    df <- .

    # MODEL 1: propensity score with covariates
    # m_prop <- multinom(rs_77924615_A_G_genotype ~ diabetes + BMI_cat + albuminuria + pgs_cross_594_umod + smoking,
    #                    data = df, trace = FALSE)
    m_prop <- multinom(rs_77924615_A_G_genotype ~ diabetes + ns(sc_BMI, df=4) + ns(sc_UACR, df=4) + ns(pgs_cross_594_umod, df=4) + smoking,
                       data = df, trace = FALSE)
    prob_mat_prop <- predict(m_prop, newdata = df, type = "prob")

    # MODEL 2: raw genotype probability (no covariates)
    m_raw <- multinom(rs_77924615_A_G_genotype ~ 1, data = df, trace = FALSE)
    prob_mat_raw <- predict(m_raw, newdata = df, type = "prob")

    # Index for observed genotype
    geno_index <- as.integer(df$rs_77924615_A_G_genotype)

    # Extract scores
    pscore <- vapply(seq_len(nrow(df)), function(i) prob_mat_prop[i, geno_index[i]], numeric(1))
    rawscore <- vapply(seq_len(nrow(df)), function(i) prob_mat_raw[i, geno_index[i]], numeric(1))

    df %>% mutate(
      propensity_score_0 = prob_mat_prop[, 1],
      propensity_score_1 = prob_mat_prop[, 2],
      propensity_score_2 = prob_mat_prop[, 3],
      propensity_score = pscore,
      raw_score        = rawscore,
      stabilized_weight = raw_score / propensity_score
    )
  }) %>%
  ungroup()

events_ps <- events_ps %>%
  mutate(
    genotype = factor(
      rs_77924615_A_G_genotype,
      levels = c(0, 1, 2),
      labels = c("0", "1", "2")
    )
  ) %>%
  mutate(from_num = gsub("State ", "", from)) %>%
  mutate(from_num = as.character(from_num))  # match ped_events$from

### descriptive plots (old) ----

# Generalized function for binary variables (1 = "yes")
binom_ci_tbl <- function(data, var, group = from) {
  out <- data %>%
    group_by({{group}}) %>%
    summarise(
      n = sum(!is.na(!!sym(var))),
      x = sum(!!sym(var) == 1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      binom = purrr::map2(x, n, ~ binom.test(.x, .y)),
      prop  = x / n,
      lower = purrr::map_dbl(binom, ~ .x$conf.int[1]),
      upper = purrr::map_dbl(binom, ~ .x$conf.int[2])
    ) %>%
    select({{group}}, prop, lower, upper) %>%
    mutate(
      prop  = percent(prop,  accuracy = 0.1),
      lower = percent(lower, accuracy = 0.1),
      upper = percent(upper, accuracy = 0.1),
      ci    = paste0("(", lower, ", ", upper, ")"),
      variable = var
    )

  out <- out %>%
    select(-c(lower, upper)) %>%
    arrange({{group}})

  return(out)
}

tab_diabetes    <- binom_ci_tbl(events_w, "diabetes")
tab_smoking     <- binom_ci_tbl(events_w, "smoking")

library(dplyr)
library(purrr)
library(scales)

cat_percent_ci_tbl <- function(data, var, group = "from") {
  vals <- na.omit(unique(data[[var]]))
  n_levels <- length(vals)
  group_syms <- syms(group)
  out <- NULL

  # if (n_levels == 2) {
  #   out <- data %>%
  #     group_by(!!!group_syms) %>%
  #     summarise(
  #       n = sum(!is.na(.data[[var]])),
  #       x = sum(.data[[var]] == vals[2], na.rm = TRUE),
  #       .groups = "drop"
  #     ) %>%
  #     mutate(
  #       binom = map2(x, n, ~ binom.test(.x, .y)),
  #       prop  = x / n,
  #       lower = map_dbl(binom, ~ .x$conf.int[1]),
  #       upper = map_dbl(binom, ~ .x$conf.int[2]),
  #       level = vals[2],
  #       variable = var,
  #       prop  = percent(prop,  accuracy = 0.1),
  #       lower = percent(lower, accuracy = 0.1),
  #       upper = percent(upper, accuracy = 0.1),
  #       ci    = paste0("(", lower, ", ", upper, ")")
  #     ) %>%
  #     select(!!!group_syms, level, prop, variable, ci)
  # } else {
    out <- map_dfr(vals, function(vv) {
      data %>%
        group_by(!!!group_syms) %>%
        summarise(
          n = sum(!is.na(.data[[var]])),
          x = sum(.data[[var]] == vv, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        mutate(
          binom = map2(x, n, ~ binom.test(.x, .y)),
          prop  = x / n,
          lower = map_dbl(binom, ~ .x$conf.int[1]),
          upper = map_dbl(binom, ~ .x$conf.int[2]),
          level = vv,
          variable = var,
          prop  = percent(prop,  accuracy = 0.1),
          lower = percent(lower, accuracy = 0.1),
          upper = percent(upper, accuracy = 0.1),
          ci    = paste0("(", lower, ", ", upper, ")")
        )
    }) %>%
      select(variable, !!!group_syms, level, prop, ci)
  # }
  out %>% arrange(across(all_of(group)), level)
}

tab_albuminuria <- cat_percent_ci_tbl(events_w, "albuminuria")
tab_bmi <- cat_percent_ci_tbl(events_w, "BMI_cat")

genotypes <- levels(events_w2$rs_77924615_A_G_genotype)

tab_diabetes_by_genotype <- purrr::map_dfr(genotypes, function(gt) {
  cat_percent_ci_tbl(
    data = dplyr::filter(events_w2, rs_77924615_A_G_genotype == gt),
    var = "diabetes",
    group = "from"
  ) %>%
    mutate(genotype = gt)
})
tab_diabetes_by_genotype %>% filter(level == 1)

tab_albuminuria_by_genotype <- purrr::map_dfr(genotypes, function(gt) {
  cat_percent_ci_tbl(
    data = dplyr::filter(events_w2, rs_77924615_A_G_genotype == gt),
    var = "albuminuria",
    group = "from"
  ) %>%
    mutate(genotype = gt)
})
tab_albuminuria_by_genotype %>% print(n=30)

tab_bmi_by_genotype <- purrr::map_dfr(genotypes, function(gt) {
  cat_percent_ci_tbl(
    data = dplyr::filter(events_w2, rs_77924615_A_G_genotype == gt),
    var = "BMI_cat",
    group = "from"
  ) %>%
    mutate(genotype = gt)
})
tab_bmi_by_genotype %>% print(n=30)

tab_smoking_by_genotype <- purrr::map_dfr(genotypes, function(gt) {
  cat_percent_ci_tbl(
    data = dplyr::filter(events_w2, rs_77924615_A_G_genotype == gt),
    var = "smoking",
    group = "from"
  ) %>%
    mutate(genotype = gt)
})
tab_smoking_by_genotype %>% filter(level == 1)

### UMOD dosages by state
hist_umod <- ggplot(events_w, aes(
    x      = rs77924615_A_G,
    weight = weight
  )) +
  geom_histogram(
    bins     = 30,
    color    = "black",
    fill     = "steelblue",
    boundary = 0
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "UMOD dosages by state",
    x     = "rs77924615_A_G",
    y     = "Relative frequency"
  )

ggsave(file.path(dir_out, "/figures/descriptives/UMOD_histogram_states_combined.png"), plot = hist_umod, width = 10, height = 18, units = "cm")

### PSGs (594) by state
hist_pgs_594_cross <- ggplot(events_w, aes(
    x      = pgs_cross_594_umod,
    weight = weight
  )) +
  geom_histogram(
    bins     = 30,
    color    = "black",
    fill     = "steelblue",
    boundary = 0
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "PGS 594 cross by state",
    x     = "PGS 594 cross",
    y     = "Relative frequency"
  )

ggsave(file.path(dir_out, "/figures/descriptives/PGS594cross_histogram_states_combined.png"), plot = hist_pgs_594_cross, width = 10, height = 18, units = "cm")

### diabetes status by state
barplot_diabetes <- ggplot(events_w, aes(
    x      = diabetes,
    weight = weight
  )) +
  geom_bar(
    color    = "black",
    fill     = "steelblue"
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Diabetes status by state",
    x     = "Diabetes status",
    y     = "Relative frequency"
  )

ggsave(file.path(dir_out, "/figures/descriptives/diabetes_barplot_states_combined.png"), plot = barplot_diabetes, width = 10, height = 18, units = "cm")

### diabetes status by state and UMOD genotype
for (genotype in c("0", "1", "2")) {
  barplot_diabetes_genotype <- ggplot(
    events_w2 %>% filter(rs_77924615_A_G_genotype == genotype),
    aes(x = diabetes, weight = weight_g)
  ) +
    geom_bar(color = "black", fill = "steelblue") +
    facet_wrap(~ from, ncol = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste("Diabetes for UMOD genotype", genotype),
      x     = "Diabetes status",
      y     = "Relative frequency"
    )

  ggsave(file.path(dir_out, paste0("/figures/descriptives/diabetes_barplot_states_combined_genotype_", genotype, ".png")), plot = barplot_diabetes_genotype, width = 10, height = 18, units = "cm")
}

### UACR levels by state
hist_uacr <- ggplot(events_w, aes(
    x      = sc_UACR,
    weight = weight
  )) +
  geom_histogram(
    bins     = 30,
    color    = "black",
    fill     = "steelblue",
    boundary = 0
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "UACR levels by state",
    x     = "UACR",
    y     = "Relative frequency"
  )
ggsave(file.path(dir_out, "/figures/descriptives/UACR_histogram_states_combined.png"), plot = hist_uacr, width = 10, height = 18, units = "cm")

# albuminuria by state
barplot_albuminuria <- ggplot(events_w, aes(
    x      = albuminuria,
    weight = weight
  )) +
  geom_bar(
    color    = "black",
    fill     = "steelblue"
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Albuminuria by state",
    x     = "Albuminuria",
    y     = "Relative frequency"
  )
ggsave(file.path(dir_out, "/figures/descriptives/albuminuria_barplot_states_combined.png"), plot = barplot_albuminuria, width = 10, height = 18, units = "cm")

# albuminuria by state and UMOD genotype
for (genotype in c("0", "1", "2")) {
  barplot_albuminuria_genotype <- ggplot(
    events_w2 %>% filter(rs_77924615_A_G_genotype == genotype),
    aes(x = albuminuria, weight = weight_g)
  ) +
    geom_bar(color = "black", fill = "steelblue") +
    facet_wrap(~ from, ncol = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste("Albuminuria for UMOD genotype", genotype),
      x     = "Albuminuria",
      y     = "Relative frequency"
    )

  ggsave(file.path(dir_out, paste0("/figures/descriptives/albuminuria_barplot_states_combined_genotype_", genotype, ".png")), plot = barplot_albuminuria_genotype, width = 10, height = 18, units = "cm")
}

# BMI by state
hist_bmi <- ggplot(events_w, aes(
    x      = sc_BMI,
    weight = weight
  )) +
  geom_histogram(
    bins     = 30,
    color    = "black",
    fill     = "steelblue",
    boundary = 0
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "BMI levels by state",
    x     = "BMI",
    y     = "Relative frequency"
  )
ggsave(file.path(dir_out, "/figures/descriptives/BMI_histogram_states_combined.png"), plot = hist_bmi, width = 10, height = 18, units = "cm")

# BMI categories by state
barplot_bmi_cat <- ggplot(events_w, aes(
    x      = BMI_cat,
    weight = weight
  )) +
  geom_bar(
    color    = "black",
    fill     = "steelblue"
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "BMI categories by state",
    x     = "BMI category",
    y     = "Relative frequency"
  )
ggsave(file.path(dir_out, "/figures/descriptives/BMI_cat_barplot_states_combined.png"), plot = barplot_bmi_cat, width = 10, height = 18, units = "cm")

# BMI categories by state and UMOD genotype
for (genotype in c("0", "1", "2")) {
  barplot_bmi_cat_genotype <- ggplot(
    events_w2 %>% filter(rs_77924615_A_G_genotype == genotype),
    aes(x = BMI_cat, weight = weight_g)
  ) +
    geom_bar(color = "black", fill = "steelblue") +
    facet_wrap(~ from, ncol = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste("BMI cat for UMOD genotype", genotype),
      x     = "BMI category",
      y     = "Relative frequency"
    )

  ggsave(file.path(dir_out, paste0("/figures/descriptives/BMI_cat_barplot_states_combined_genotype_", genotype, ".png")), plot = barplot_bmi_cat_genotype, width = 10, height = 18, units = "cm")
}

# Smoking by state
barplot_smoking <- ggplot(events_w, aes(
    x      = smoking,
    weight = weight
  )) +
  geom_bar(
    color    = "black",
    fill     = "steelblue"
  ) +
  facet_wrap(~ from, ncol = 1) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Smoking status by state",
    x     = "Smoking status",
    y     = "Relative frequency"
  )
ggsave(file.path(dir_out, "/figures/descriptives/smoking_barplot_states_combined.png"), plot = barplot_smoking, width = 10, height = 18, units = "cm")

# Smoking by state and UMOD genotype
for (genotype in c("0", "1", "2")) {
  barplot_smoking_genotype <- ggplot(
    events_w2 %>% filter(rs_77924615_A_G_genotype == genotype),
    aes(x = smoking, weight = weight_g)
  ) +
    geom_bar(color = "black", fill = "steelblue") +
    facet_wrap(~ from, ncol = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(
      title = paste("Smoking for UMOD genotype", genotype),
      x     = "Smoking status",
      y     = "Relative frequency"
    )

  ggsave(file.path(dir_out, paste0("/figures/descriptives/smoking_barplot_states_combined_genotype_", genotype, ".png")), plot = barplot_smoking_genotype, width = 10, height = 18, units = "cm")
}

## propensity scores ----

### scatter plots ----
# Define transition pairs
transitions <- list("0_1" = c(0, 1), "1_2" = c(1, 2))

# Define per‐geno axis ranges
ps_ranges <- list(
  "0" = c(0.0, 0.2),
  "1" = c(0.1, 0.5),
  "2" = c(0.4, 0.8)
)

# Loop over transitions and genotype levels
for (trans_name in names(transitions)) {
  states <- transitions[[trans_name]]
  from_k <- states[1]
  to_k   <- states[2]

  # Filter and pivot data
  e_pair <- events_ps %>%
    filter(from_num %in% states) %>%
    group_by(id) %>%
    filter(n_distinct(from_num) == 2) %>%
    ungroup() %>%
    select(id, from_num, starts_with("propensity_score_")) %>%
    pivot_wider(
      names_from  = from_num,
      values_from = starts_with("propensity_score_"),
      names_sep   = "_from_"
    )

  # Loop over genotype levels
  for (geno in 0:2) {
    col_x <- paste0("propensity_score_", geno, "_from_", from_k)
    col_y <- paste0("propensity_score_", geno, "_from_", to_k)

    # Safety check
    if (!all(c(col_x, col_y) %in% colnames(e_pair))) next

    # grab the appropriate axis limits
    lims <- ps_ranges[[as.character(geno)]]

    # Create plot
    p <- ggplot(e_pair, aes_string(x = col_x, y = col_y)) +
      # much smaller pure‐black dots, no outline/shadow:
      geom_point(size = 0.7, shape = 16, color = "black", alpha = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste0("PS_", geno, " ", from_k, "→", to_k),
        x     = paste0("Propensity Score ", geno, " (from = ", from_k, ")"),
        y     = paste0("Propensity Score ", geno, " (from = ", to_k, ")")
      ) +
      coord_fixed() +
      scale_x_continuous(
        breaks = seq(lims[1], lims[2], by = 0.1),
        limits = lims
      ) +
      scale_y_continuous(
        breaks = seq(lims[1], lims[2], by = 0.1),
        limits = lims
      )

    # Save
    filename <- paste0("ps_scatter_", trans_name, "_geno", geno, ".png")
    ggsave(
      file.path(dir_out, "figures/descriptives/ps", filename),
      plot = p,
      width = 8, height = 8, units = "cm"
    )
  }
}

### models ----

ped_events_with_weights <- ped_events %>%
  left_join(
    events_ps %>% select(id, from_num, stabilized_weight),
    by = c("id", "from" = "from_num")
  )

# 1. Define transitions/from codes
transitions <- tibble(
  transition = c("0->1", "1->2", "2->3"),
  from_value = c("0",    "1",    "2")
)

# 2. Build base and adjustment terms
base_terms <- "s(tend, by = transition) + sex*transition + rs77924615_A_G*transition"
adj_terms  <- paste0(
  "pgs_cross_594_umod*transition + diabetes*transition + smoking*transition + ",
  "s(sc_BMI, by = transition) + s(sc_UACR, by = transition) + s(eGFRcrea, by = transition)"
)

# 3. Helper to generate per-transition formulas
make_formulas <- function(tr) {
  extras <- switch(
    tr,
    "0->1" = "",
    "1->2" = " + s(age_onset, by = transition_after_onset_strat)",
    "2->3" = paste0(
      " + s(age_onset, by = transition_after_onset_strat)",
      " + s(age_progression, by = transition_after_progression)"
    )
  )
  list(
    unadj = as.formula(paste0("ped_status ~ ", base_terms, extras)),
    adj_only_eGFR = as.formula(paste0("ped_status ~ ", base_terms, extras, " + s(eGFRcrea, by = transition)")),
    adj   = as.formula(paste0("ped_status ~ ", base_terms, " + ", adj_terms, extras))
  )
}

# pre-build all formulas
formula_list <- set_names(
  map(transitions$transition, make_formulas),
  transitions$transition
)

# 1) Build your transition grid as before
grid <- expand.grid(
  transition   = c("0->1", "1->2", "2->3"),
  formula_type = c("unadj", "adj_only_eGFR", "adj"),
  weight_type  = c("unw",   "w"),
  stringsAsFactors = FALSE
)

# 2) Augment it with the other info you need
grid <- grid %>%
  mutate(
    from_value = case_when(
      transition == "0->1" ~ "0",
      transition == "1->2" ~ "1",
      transition == "2->3" ~ "2"
    ),
    variant = paste0(formula_type, "_", weight_type),
    weight_var = if_else(weight_type == "w", "stabilized_weight", NA_character_)
  )

# 3) Set up your parallel backend
cl <- makeCluster(18)
registerDoParallel(cl)

# 4) Run all 12 fits in one go, *passing* each column of `grid` in as its own argument
results <- foreach(
  transition   = grid$transition,
  from_value   = grid$from_value,
  formula_type = grid$formula_type,
  weight_type  = grid$weight_type,
  variant      = grid$variant,
  weight_var   = grid$weight_var,
  .combine     = bind_rows,
  .packages    = c("dplyr","mgcv")
) %dopar% {

  # a) subset the data
  dat <- ped_events_with_weights %>% filter(from == from_value)

  # b) build the right formula on the fly
  base_terms <-
    "s(tend, by = transition) + sex*transition + rs77924615_A_G*transition"
  adj_terms  <- paste0(
    "pgs_cross_594_umod*transition + diabetes*transition + smoking*transition + ",
    "s(sc_BMI, by = transition) + s(sc_UACR, by = transition) + s(eGFRcrea, by = transition)"
  )
  extras <- switch(
    transition,
    "0->1" = "",
    "1->2" = " + s(age_onset, by = transition_after_onset_strat)",
    "2->3" = paste0(
      " + s(age_onset, by = transition_after_onset_strat)",
      " + s(age_progression, by = transition_after_progression)"
    )
  )

  f <- if (formula_type == "unadj") {
    as.formula(paste0("ped_status ~ ", base_terms, extras))
  } else {
    as.formula(paste0("ped_status ~ ", base_terms, " + ", adj_terms, extras))
  }

  # c) pick up the weights vector if requested
  wvec <- if (is.na(weight_var)) NULL else dat[[weight_var]]

  # d) fit the model
  m <- bam(
    formula  = f,
    data     = dat,
    family   = poisson(),
    offset   = dat$offset,
    method   = "fREML",
    discrete = TRUE,
    weights  = wvec
  )

  # e) extract the rs77924615_A_G row
  pt <- summary(m)$p.table["rs77924615_A_G", ]

  # f) return a one‐row tibble
  tibble(
    transition = transition,
    variant    = variant,
    coef       = pt["Estimate"],
    se         = pt["Std. Error"],
    p          = pt["Pr(>|z|)"]
  )
}

# 5) clean up
stopCluster(cl)

out <- results %>%
  mutate(
    coef = round(coef, 2),
    se = round(se, 2),
    p = formatC(p, format = "e", digits = 2)
  )

write.csv2(out, file = file.path(dir_out, "weighting_and_adjusting_models_results.csv"))

View(results %>% arrange(variant, transition))

### histogram ----

# hist_ps <- ggplot(events_ps, aes(
#     x     = propensity_score,
#     color = genotype,
#     fill  = genotype,
#     group = genotype
#   )) +
#   geom_density(
#     aes(y = after_stat(density)),
#     alpha = 0.3,
#     linewidth  = 1
#   ) +
#   facet_wrap(~ from, ncol = 1) +
#   scale_x_continuous(
#     name   = "Propensity score",
#     limits = c(0, 1)
#   ) +
#   scale_y_continuous(
#     name   = "Density"
#   ) +
#   scale_color_brewer(
#     palette = "Set1",
#     name    = "Genotype"
#   ) +
#   scale_fill_brewer(
#     palette = "Set1",
#     name    = "Genotype"
#   ) +
#   labs(
#     title = "In-Sample Propensity Score Densities\nby UMOD Genotype and State"
#   ) +
#   theme(
#     strip.text = element_text(face = "bold"),
#     legend.position = "bottom"
#   )

# ggsave(file.path(dir_out, "/figures/descriptives/UMOD_ps_histogram_states_combined.png"), plot = hist_ps, width = 10, height = 18, units = "cm")

## bar plots by risk factor and genotype ----

risk_factors <- c("diabetes", "albuminuria", "BMI_cat", "smoking")

### unweighted

# for (risk_factor in risk_factors) {

#   p <- ggplot(
#     events_ps %>% filter(!is.na(.data[[risk_factor]])),
#     aes_string(x = risk_factor, weight = "weight_g")
#   ) +
#     geom_bar(color = "black", fill = "steelblue", width = 0.7) +
#     # rows = from, cols = genotype
#     facet_grid(from ~ rs_77924615_A_G_genotype) +
#     scale_y_continuous(labels = percent_format(accuracy = 1)) +
#     labs(
#       title = paste(risk_factor, " for UMOD genotype by state"),
#       x     = risk_factor,
#       y     = "Relative frequency"
#     ) +
#     theme(
#       strip.background = element_rect(fill = "grey95"),
#       strip.text      = element_text(face = "bold")
#     )

#   out_file <- file.path(
#     dir_out,
#     "figures/descriptives/barplots",
#     paste0("unweighted_", risk_factor, ".png")
#   )
#   ggsave(out_file, plot = p, width = 18, height = 18, units = "cm")
# }

for (risk_factor in risk_factors) {

  p <- ggplot(
    events_ps %>% filter(!is.na(.data[[risk_factor]])),
    aes_string(x = risk_factor)
  ) +
    geom_bar(
      aes(
        y     = after_stat(prop),
        group = 1
      ),
      stat    = "count",
      color   = "black",
      fill    = "steelblue",
      width   = 0.7
    ) +
    facet_grid(from ~ rs_77924615_A_G_genotype) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(
      title = paste(risk_factor, "for UMOD genotype by state", "(unweighted)"),
      x     = risk_factor,
      y     = "Relative frequency"
    ) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text       = element_text(face = "bold")
    )

  ggsave(
    file.path(
      dir_out,
      "figures/descriptives/barplots",
      paste0("unweighted_", risk_factor, ".png")
    ),
    plot  = p,
    width = 18, height = 18, units = "cm"
  )
}

## weighted by stabilized_weight
for (risk_factor in risk_factors) {

  p <- ggplot(
    events_ps %>% filter(!is.na(.data[[risk_factor]])),
    aes_string(x = risk_factor, weight = "stabilized_weight")
  ) +
    geom_bar(
      aes(
        y     = after_stat(prop),
        group = 1
      ),
      stat    = "count",
      color   = "black",
      fill    = "steelblue",
      width   = 0.7
    ) +
    facet_grid(from ~ rs_77924615_A_G_genotype) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(
      title = paste(risk_factor, "for UMOD genotype by state", "(weighted)"),
      x     = risk_factor,
      y     = "Relative frequency"
    ) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text       = element_text(face = "bold")
    )

  ggsave(
    file.path(
      dir_out,
      "figures/descriptives/barplots",
      paste0("weighted_", risk_factor, ".png")
    ),
    plot  = p,
    width = 18, height = 18, units = "cm"
  )
}

### SMD ----

compute_smd <- function(df,
                        risk_factor,
                        state,
                        genotype_levels = c(0, 1),
                        weights = NULL) {
  # 1) Subset and drop NAs
  sub <- df %>%
    filter(from == state,
           rs_77924615_A_G_genotype %in% genotype_levels,
           !is.na(.data[[risk_factor]]),
           !is.na(rs_77924615_A_G_genotype)) %>%
    { if (!is.null(weights)) filter(., !is.na(.data[[weights]])) else . }

  # 2) Extract vector and weights
  vec  <- sub[[risk_factor]]
  wvec <- if (is.null(weights)) rep(1, nrow(sub)) else sub[[weights]]

  # 3) Factorify genotype
  g    <- factor(sub$rs_77924615_A_G_genotype, levels = genotype_levels)
  n1   <- sum(wvec[g==genotype_levels[1]])
  n2   <- sum(wvec[g==genotype_levels[2]])

  # 4) Branch on type of vec
  if (is.numeric(vec) || is.integer(vec)) {
    # — continuous case (same as before) —
    x    <- as.numeric(vec)
    summs <- tibble(x, g, w = wvec) %>%
      group_by(g) %>%
      summarise(
        m  = sum(w * x) / sum(w),
        sd = sqrt(sum(w * (x - m)^2) / sum(w)),
        .groups = "drop"
      ) %>%
      arrange(g)

    m1 <- summs$m[1]; m2 <- summs$m[2]
    sd1 <- summs$sd[1]; sd2 <- summs$sd[2]
    sd_pooled <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
    return((m1 - m2) / sd_pooled)

  } else {
    # treat as categorical
    fac <- factor(vec)
    levs <- levels(fac)

    if (length(levs) == 2) {
      # — binary factor case —
      # convert to 0/1 indicator where second level == 1
      ind <- as.integer(fac == levs[2])
      p1  <- weighted.mean(ind[g==genotype_levels[1]], wvec[g==genotype_levels[1]])
      p2  <- weighted.mean(ind[g==genotype_levels[2]], wvec[g==genotype_levels[2]])
      # pooled prop
      p   <- (n1 * p1 + n2 * p2) / (n1 + n2)
      return((p1 - p2) / sqrt(p * (1 - p)))

    } else {
      stop("`", risk_factor, "` has ", length(levs),
           " levels; only continuous or binary supported.")
    }
  }
}

risk_factors <- c("pgs_cross_594_umod", "diabetes", "sc_BMI", "sc_UACR", "smoking", "eGFRcrea")
states <- c("State 0", "State 1", "State 2")
genotype_combos <- list(
  c(0, 1), # UMOD genotype 0 and 1
  c(0, 2), # UMOD genotype 0 and 2
  c(1, 2)  # UMOD genotype 1 and 2
)
combo_tbl <- tibble(genotype_levels = genotype_combos)
weights_list <- list(NULL, "stabilized_weight")

smd_results <- expand_grid(
  risk_factor      = risk_factors,
  state            = states,
  genotype_levels  = genotype_combos,
  weights          = weights_list
) %>%
  # create clean labels
  mutate(
    genotype_levels = map_chr(genotype_levels, ~ paste0(.x, collapse = "_")),
    weights         = if_else(map_lgl(weights, is.null),
                              "none",
                              "stabilized_weight")
  ) %>%
  rowwise() %>%
  mutate(
    SMD = compute_smd(
      df               = events_ps,
      risk_factor      = risk_factor,
      state            = state,
      genotype_levels  = as.integer(strsplit(genotype_levels, "_")[[1]]),
      weights          = if (weights == "none") NULL else "stabilized_weight"
    )
  ) %>%
  ungroup() %>%
  select(risk_factor, state, genotype_levels, weights, SMD)

#### table ----
smd_results

#### love plot ----
plot_love <- function(smd_df, genotype_lvl, st) {
  # 1) Filter to the combo
  df <- smd_df %>%
    filter(
      genotype_levels == genotype_lvl,
      state            == st
    ) %>%
    mutate(
      abs_SMD = abs(SMD),
      # Preserve the exact risk‐factor order as they appear in smd_df
      risk_factor = factor(
        risk_factor,
        levels = unique(risk_factor)
      )
    ) %>%
    # 2) Sort by weights, then by risk_factor so the path goes straight down
    arrange(weights, risk_factor)

  ggplot(df, aes(
      x     = abs_SMD,
      y     = risk_factor,
      group = weights,    # one ribbon per weight scheme
      color = weights
    )) +
    # vertical connectors
    geom_path(size = 1) +
    # points
    geom_point(size = 3) +
    scale_color_manual(
      values = c(none = "red", stabilized_weight = "blue"),
      labels = c(none = "Unweighted", stabilized_weight = "Weighted")
    ) +
    # 0.1 threshold
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +
    labs(
      x     = "Absolute Standardized Mean Difference",
      y     = NULL,
      color = "Weighting",
      title = paste0(
        "Love plot - Genotype: ", genotype_lvl,
        " | State: ", st
      )
    ) +
    theme(
      panel.background    = element_blank(),
      panel.grid.major.x  = element_line(color = "grey90"),
      panel.grid.major.y  = element_blank(),
      panel.grid.minor    = element_blank(),
      axis.text.y         = element_text(size = 12),
      legend.position     = "bottom"
    )
}

genotype_levels_list <- c("0_1", "1_2", "0_2")
states               <- c("State 0", "State 1", "State 2")

# for (gen in genotype_levels_list) {
#   for (st in states) {
#     # make the plot
#     p <- plot_love(smd_results, gen, st)

#     # build a filename like "1_2State0.png"
#     fname <- paste0(gen, gsub(" ", "", st), ".png")

#     # save it
#     ggsave(
#       filename = file.path(dir_out, "figures", "descriptives", "loveplots", fname),
#       plot     = p,
#       width    = 18,
#       height   = 18,
#       units    = "cm"
#     )
#   }
# }

library(patchwork)
for (st in states) {
  # build the three plots, adding a subtitle for clarity
  plots <- lapply(genotype_levels_list, function(gen) {
    plot_love(smd_results, gen, st) +
      labs(subtitle = paste0("Genotype: ", gen))
  })

  # stitch them horizontally
  combined <- plots[[1]] | plots[[2]] | plots[[3]]

  # filename based only on state
  fname <- paste0(gsub(" ", "", st), ".png")

  ggsave(
    filename = file.path(dir_out, "figures", "descriptives", "loveplots", fname),
    plot     = combined,
    width    = 54,    # 3 × 18 cm
    height   = 18,
    units    = "cm"
  )
}

# Aalen-Johansen estimator ----
tra.idm <- matrix(FALSE, 5, 5, dimnames = list(c(0, 1, 2, 3, 4), c(0, 1, 2, 3, 4)))
tra.idm[1, ] <- c(F, T, F, F, T)
tra.idm[2, ] <- c(F, F, T, F, T)
tra.idm[3, ] <- c(F, F, F, T, T)
tra.idm[4, ] <- c(F, F, F, F, F)
tra.idm[5, ] <- c(F, F, F, F, F)

events_aj_age <- events %>%
    filter(agestart >= min_age) %>%
    rename(entry = agestart, exit = agestop) %>%
    select(id, from, to, entry, exit)
state.names <- c("0", "1", "2", "3", "4")
transitions <- c("0 1", "1 2", "2 3", "0 4", "1 4", "2 4")
etm.idm <- etm(data = events_aj_age, state.names = state.names, tra = tra.idm, cens.name = "cens", s = 35)

### plot
png(file.path(dir_out, "/figures/aj_age.png"))
par(mfrow=c(2,3))
for(i in 1:length(transitions)){
  plot(etm.idm, tr.choice = transitions[i], conf.int = TRUE,
      lwd = 2, legend = FALSE, ylim = c(0, 1),
      xlim = c(35, 80), xlab = "Years",
      ci.fun = "cloglog", main = transitions[i])
}
dev.off()

## replicate using PAMs
events_ct <- events %>%
    rename(tstart = agestart, tstop = agestop) %>%
    add_counterfactual_transitions()

ped_aj <- as_ped(
  data         = events_ct,
  formula      = Surv(tstart, tstop, status) ~ .,
#   cut          = seq(19, 80, 1),
  transition   = "transition",
  id           = "id",
  censor_code  = 0,
  timescale    = "calendar") %>%
  filter(tstart >= min_age) # important!!!

formula_aj <- ped_status ~ s(tend, by = transition, k=10) + transition

pam_aj <- mgcv::bam(
  formula = formula_aj,
  data = ped_aj,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_aj)

newData_aj <- make_newdata(ped_aj, tend = sort(unique(tend)), transition=unique(transition)) %>%
  select(tstart, tend, intlen, interval, id, offset, ped_status, transition) %>%
  group_by(transition) %>%
  add_hazard(pam_aj) %>%
  add_cumu_hazard(pam_aj) %>%
  add_trans_prob(pam_aj) %>%
  ungroup()

survCurves_aj <- ggplot(newData_aj, aes(x=tend, y=trans_prob)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/pam_aj_survCurves.png"), plot = survCurves_aj, width = 10, height = 10, units = "cm")

# ### table
# # summary(etm.idm)$"0 1"[, c("P", "lower", "upper")]

# ### landmark method for probability plots (given from status at time s, probability of being in state to at time t, s <= t)
# time.points <- TBD
# landmark.etm <- lapply(time.points, function(start.time) {
#     etm(events_aj, c("0", "1", "2"), tra.idm, "cens", start.time)
#     })

# #### plot landmark
# TBD

# PEM & PAM ----

trans_2d    <- c("0->1", "0->4")
trans_3d    <- c("1->2", "2->3", "1->4", "2->4")
transitions <- c("0->1", "0->4", "1->2", "1->4", "2->3", "2->4")
age_onset_slices  <- c(40, 50, 60, 70)
prog_after_onset_slices <- c(2, 5, 10, 15)

formulas <- list(
  pam_stss = "ped_status ~
    s(tend, by = transition) +
    s(age_onset, by = transition_after_onset_strat) +
    s(age_progression, by = transition_after_progression) +
    sex*transition",
  pam_mts = "ped_status ~
    s(tend, by = transition_to_death) +
    s(tend_onset, by = transition_after_onset) +
    s(age_onset, by = transition_after_onset) +
    s(tend_progression, by = transition_after_progression) +
    s(age_progression, by = transition_after_progression) +
    sex*transition"
)

num_cores   <- length(formulas)

registerDoParallel(cores = num_cores)
foreach(model = names(formulas), .packages = c("ggplot2", "pammtools")) %dopar% {

  fig_dir <- file.path(dir_out, "figures", model)
  if (!dir.exists(fig_dir)) {
    dir.create(fig_dir, recursive = TRUE)
  }

  mod_formula <- as.formula(formulas[[model]])

  mod <- mgcv::bam(
    formula = mod_formula,
    data = ped_events,
    family=poisson(),
    offset=offset,
    method="fREML",
    discrete = TRUE)

  # ped_new <- ped_events %>% make_newped(mod)
  # saveRDS(ped_new, file.path(dir_data, sprintf("ped_%s.rds", model)))
  ped_new <- readRDS(file.path(dir_data, sprintf("ped_%s.rds", model)))

  for(trans in trans_2d) {
    trans_print <- gsub("->", "", trans)
    plots <- create_2d_plots(ped_new, model, trans)
    saveRDS(plots, file.path(dir_data, paste0(model, "_plots_", trans_print, ".rds")))
  }

  for(trans in trans_3d) {
    trans_print <- gsub("->", "", trans)
    plots_age <- create_3d_plots(ped_new, model, trans, "age", age_onset_slices, prog_after_onset_slices)
    saveRDS(plots_age, file.path(dir_data, paste0(model, "_plots_age_", trans_print, ".rds")))
    plots_time <- create_3d_plots(ped_new, model, trans, "time", age_onset_slices, prog_after_onset_slices)
    saveRDS(plots_time, file.path(dir_data, paste0(model, "_plots_time_", trans_print, ".rds")))
  }

  # lapply(trans_2d, function(trans) {create_2d_plots(ped_new, model, trans)})
  # lapply(trans_3d, function(trans) {create_3d_plots(ped_new, model, trans, "age", age_onset_slices, prog_after_onset_slices)})
  # lapply(trans_3d, function(trans) {create_3d_plots(ped_new, model, trans, "time", age_onset_slices, prog_after_onset_slices)})

}
stopImplicitCluster()

## null models ----
### pam stss ----
formula_stss <- "ped_status ~
  s(tend, by = transition) +
  s(age_onset, by = transition_after_onset_strat) +
  s(age_progression, by = transition_after_progression) +
  sex*transition"

pam_stss <- mgcv::bam(
  formula = as.formula(formula_stss),
  data = ped_events,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

summary(pam_stss)

### pam mts ----
formula_mts <- "ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  sex*transition"

pam_mts <- mgcv::bam(
  formula = as.formula(formula_mts),
  data = ped_events,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

summary(pam_mts)

AIC(pam_stss, pam_mts)
anova(pam_stss, pam_mts)

## genetic models ----

### pam stss ----
formula_stss_genetics_naive <- "ped_status ~
  s(tend, by = transition) +
  s(age_onset, by = transition_after_onset_strat) +
  s(age_progression, by = transition_after_progression) +
  sex*transition +
  rs77924615_A_G*transition"

pam_stss_genetics_naive <- mgcv::bam(
  formula = as.formula(formula_stss_genetics_naive),
  data = ped_events,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

summary(pam_stss_genetics_naive)

### pam mts ----
formula_mts_genetics_naive <- "ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  sex*transition +
  rs77924615_A_G*transition"

pam_mts_genetics_naive <- mgcv::bam(
  formula = as.formula(formula_mts_genetics_naive),
  data = ped_events,
  family = poisson(),
  offset = offset,
  method = "fREML",
  discrete = TRUE)

summary(pam_mts_genetics_naive)

















## null model ----

### fit model
formula_null <- ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  transition * sex

pam_null <- mgcv::bam(
  formula = formula_null,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_null)
AIC(pam_null, pam_null_strat_adj)


## null model (reduced after visual inspection) ----

### fit model
formula_null_red <- ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  age_onset*transition_after_onset +
  s(tend_progression, by = transition_after_progression) +
  transition + sex*transition

# formula_null <- ped_status ~ s(tend, by=rs77924615_A_G_genotype) + transition + sex + s(tend_onset, by=rs77924615_A_G_genotype) + s(tend_progression) + s(eGFRcrea, by=transition) + age_onset + age_progression + rs77924615_A_G*transition
# formula_null <- ped_status ~ s(tend, bs='fs', group=transition) + transition
# formula_null <- ped_status ~ s(tend, by=transition) + transition

pam_null_red <- mgcv::bam(
  formula = formula_null_red,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_null_red)

## null model reduced ----

model <- names(models)[2]

mod_formula <- as.formula(models[[model]])

mod <- mgcv::bam(
  formula = as.formula(models[[model]]),
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(mod)

ped_new <- ped_events %>% make_newped(mod)
saveRDS(ped_new, file.path(dir_data, sprintf("ped_%s.rds", model)))
ped_new <- readRDS(file.path(dir_data, sprintf("ped_%s.rds", model)))

formula_from0 = "ped_status ~ s(tend, by = transition) + transition"

mod_from0 <- mgcv::bam(
  formula = as.formula(formula_from0),
  data = ped_events %>% filter(from == 0),
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(mod_from0)

ci <- TRUE
trans_from0 <- c("0->1", "0->4")
ped_new_from0 <- ped_events %>%
    make_newdata(
        transition              = trans_from0,
        tend                    = sort(unique(tend))
    ) %>%
    mutate(
        from = as.integer(sub("([0-9])->.*","\\1", transition)),
        to   = as.integer(sub(".*->([0-9])","\\1", transition)),
        intlen = tend - tstart
    ) %>%
    add_transVars() %>%
    select(one_of(c("tstart", "tend", "from", "to", "transition", "intlen"))) %>%
    group_by(transition) %>%
    arrange(transition, tend) %>%
    add_hazard(mod_from0, type = "link", ci = ci) %>%
    rename(loghazard = hazard, loghazard_lower = ci_lower, loghazard_upper = ci_upper) %>%
    add_hazard(mod_from0, ci = ci) %>%
    rename(hazard_lower = ci_lower, hazard_upper = ci_upper) %>%
    add_cumu_hazard(mod_from0, ci = ci) %>%
    add_trans_prob(mod_from0, ci = ci) %>%
    ungroup() %>%
    mutate(transition = factor(transition, levels = c("0->1", "0->4")))

lapply(trans_from0, function(trans) {create_2d_plots(ped_new_from0, "pam_null_from0", trans)})

### likelihood ratio test reduced versus full ----
anova(pam_null, pam_null_reduced, test = "LRT")

png(file.path(dir_out, "/figures/pam_null_hazards.png"))
plot(pam_null, xlim=c(40,80), ylim=c(-5,5), page =1)
dev.off()

for(i in c(1:2, 4:5, 7:8, 10:11, 13:14)){
  png(file.path(dir_out, paste0("/figures/pam_null_hazards_", i, ".png")))
  plot.gam(pam_null, select = i)
  dev.off()
}

## null model (age stratified by transition, debugging) ----
events_ct <- events %>%
    rename(tstart = agestart, tstop = agestop) %>%
    add_counterfactual_transitions()

ped_aj <- as_ped(
  data         = events_ct,
  formula      = Surv(tstart, tstop, status) ~ .,
#   cut          = seq(19, 80, 1),
  transition   = "transition",
  id           = "id",
  censor_code  = 0,
  timescale    = "calendar") %>%
  filter(tstart >= min_age)

ped_aj_2 <- ped_aj %>%
  filter(from == 2) %>%
  filter(tstart >= 45)

ped_aj_2 %>%
  group_by(id, transition) %>%
  slice(n()) %>%
  group_by(transition) %>%
  summarise(mean(ped_status))

a <- ped_aj_2 %>% group_by(to, tend) %>% summarise(mean(ped_status))
print(a, n=50)
formula_aj <- ped_status ~ s(tend, by = transition) + transition

formula_aj <- ped_status ~ interval*transition

pam_aj_2 <- mgcv::gam(
  formula = formula_aj,
  data = ped_aj_2,
  family=poisson(),
  offset=offset
  # method="fREML",
  # discrete = TRUE
  )

newData_aj_2 <- make_newdata(ped_aj_2, tend = sort(unique(tend)), transition=unique(transition)) %>%
  select(tstart, tend, intlen, interval, id, offset, ped_status, transition) %>%
  arrange(tend, transition) %>%
  group_by(transition) %>%
  add_hazard(pam_aj_2) %>%
  add_cumu_hazard(pam_aj_2) %>%
  add_trans_prob(pam_aj_2) %>%
  ungroup() %>%
  mutate(transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")))

survCurves_aj_2 <- ggplot(newData_aj_2, aes(x=tend, y=trans_prob)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/debug/pam_aj_2_survCurves.png"), plot = survCurves_aj_2, width = 10, height = 10, units = "cm")

cumuHazard_aj_2 <- ggplot(newData_aj_2, aes(x=tend, y=cumu_hazard)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80))
  # ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/debug/pam_aj_2_cumuHazard.png"), plot = cumuHazard_aj_2, width = 10, height = 10, units = "cm")

hazard_aj_2 <- ggplot(newData_aj_2, aes(x=tend, y=hazard)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80))
  # ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/debug/pam_aj_2_hazard.png"), plot = hazard_aj_2, width = 10, height = 10, units = "cm")

library(tidyr)
all_cause_df = newData_aj_2 %>%
  select(tend, hazard, transition) %>%
  pivot_wider(
    names_from = transition,
    values_from = hazard
  ) %>%
  rename(eskd = "2->3", death = "2->4") %>%
  mutate(all_haz = eskd + death) %>%
  mutate(S = exp(-cumsum(all_haz)))

newData_aj_3 <- newData_aj_2 %>% left_join(all_cause_df, by = "tend") %>%
  group_by(transition) %>%
    mutate(cif = cumsum(hazard * S))

View(newData_aj_3 %>% select(tend, transition, hazard, all_haz, S, trans_prob, cif))

ev <- events %>%
  filter(from == 2)

library(survival)
km <- survfit(Surv( agestop, status) ~ 1, data = ev)

png(file.path(dir_out, "/figures/debug/km.png"))
plot(km, mark.time = TRUE)
dev.off()

# pam_aj <- mgcv::gam(
#   formula = formula_aj,
#   data = ped_events,
#   family=poisson(),
#   offset=offset)

pam_aj <- mgcv::bam(
  formula = formula_aj,
  data = ped_aj,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_aj)

newData_aj <- make_newdata(ped_aj, tend = sort(unique(tend)), transition=unique(transition)) %>%
  select(tstart, tend, intlen, interval, id, offset, ped_status, transition) %>%
  group_by(transition) %>%
  add_hazard(pam_aj) %>%
  add_cumu_hazard(pam_aj) %>%
  add_trans_prob(pam_aj) %>%
  ungroup() %>%
  mutate(transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")))

survCurves_aj <- ggplot(newData_aj, aes(x=tend, y=trans_prob)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/debug/pam_aj_survCurves.png"), plot = survCurves_aj, width = 10, height = 10, units = "cm")

## null model (age stratified by transition) ----

formula_null_strat <- ped_status ~ s(tend, by = transition) + transition * sex

pam_null_strat <- mgcv::bam(
  formula = formula_null_strat,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_null_strat)

newData_null_strat <- make_newdata(ped_events, tend = sort(unique(tend)), transition=unique(transition)) %>%
  group_by(transition) %>%
  add_hazard(pam_null_strat) %>%
  add_cumu_hazard(pam_null_strat) %>%
  add_trans_prob(pam_null_strat) %>%
  ungroup() %>%
  mutate(transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")))

survCurves_null_strat <- ggplot(newData_null_strat, aes(x=tend, y=trans_prob)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/pam_null_strat_survCurves.png"), plot = survCurves_null_strat, width = 10, height = 10, units = "cm")

png(file.path(dir_out, "/figures/pam_null_strat_hazards_1.png"))
plot.gam(pam_null_strat, select = 1, xlim=c(35,80), ylim=c(-3,5))
dev.off()

png(file.path(dir_out, "/figures/pam_null_strat_hazards_2.png"))
plot.gam(pam_null_strat, select = 2)
dev.off()

png(file.path(dir_out, "/figures/pam_null_strat_hazards_3.png"))
plot.gam(pam_null_strat, select = 3)
dev.off()

png(file.path(dir_out, "/figures/pam_null_strat_hazards_4.png"))
plot.gam(pam_null_strat, select = 4)
dev.off()

png(file.path(dir_out, "/figures/pam_null_strat_hazards_5.png"))
plot.gam(pam_null_strat, select = 5)
dev.off()

png(file.path(dir_out, "/figures/pam_null_strat_hazards_6.png"))
plot.gam(pam_null_strat, select = 6)
dev.off()

## null model (age stratified by transition, with left-truncation adjustment) ----

formula_null_strat_adj <- ped_status ~
  s(tend, by = transition) +
  s(age_onset, by = transition_after_onset_strat) +
  s(age_progression, by = transition_after_progression) +
  transition * sex

pam_null_strat_adj <- mgcv::bam(
  formula = formula_null_strat_adj,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_null_strat_adj)

# Define the indices and corresponding file names for the plots
plot_indices <- list(
  hazards = 1:6,
  ageOnset = c(9, 10, 11, 12),
  ageProgression = c(17, 18)
)

# Loop through each category and generate the plots
for (category in names(plot_indices)) {
  for (i in plot_indices[[category]]) {
    file_name <- sprintf("/figures/pam_null_strat_adj_%s_%d.png", category, i)
    png(file.path(dir_out, file_name))
    plot.gam(pam_null_strat_adj, select = i, xlim = c(35, 80), ylim = c(-3, 5))
    dev.off()
  }
}

newData_null_strat_adj <- make_newdata(ped_events, tend = sort(unique(tend)), transition=unique(transition)) %>%
  group_by(transition) %>%
  add_hazard(pam_null_strat_adj) %>%
  add_cumu_hazard(pam_null_strat_adj) %>%
  add_trans_prob(pam_null_strat_adj) %>%
  ungroup() %>%
  mutate(transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")))

survCurves_null_strat <- ggplot(newData_null_strat, aes(x=tend, y=trans_prob)) +
  geom_line() +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1))
ggsave(file.path(dir_out, "/figures/pam_null_strat_survCurves.png"), plot = survCurves_null_strat, width = 10, height = 10, units = "cm")

## null model (age stratified by transition, with left-truncation adjustment and tensor for age by time-since-onset) ----

### fit model
formula_tensor <- ped_status ~ ti(tend, by = transition) + ti(age_onset, by = transition_after_onset_strat) + ti(tend, age_onset, by = transition_after_onset_strat) + transition
formula_tensor <- ped_status ~ te(tend, age_onset, by = transition_after_onset_strat) + transition
formula_tensor <- ped_status ~ s(tend, by = transition) + s(age_onset, by = transition) + transition

## genetics model (pgs) ----

### fit model
formula_pgs <- ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  transition +
  sex*transition +
  pgs_decline*transition

pam_pgs <- mgcv::bam(
  formula = formula_pgs,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_pgs)

### plot survival curves
newData_pgs <- make_newdata(ped_events, tend = sort(unique(tend)), transition=unique(transition), pgs_decline = c(0, 12, 24)) %>%
  group_by(transition, pgs_decline) %>%
  add_cumu_hazard(pam_pgs) %>%
  add_trans_prob(pam_pgs)

survCurves_pgs <- ggplot(newData_pgs, aes(x=tend, y=trans_prob)) +
  geom_line(aes(col=as.factor(pgs_decline))) +
  # geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2) +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(file.path(dir_out, "/figures/pam_pgs_survCurves.png"), plot = survCurves_pgs, width = 10, height = 10, units = "cm")

## genetics model age-stratified (pgs) ----

### fit model
formula_pgs_strat <- ped_status ~
  s(tend, by = transition) +
  transition +
  sex*transition +
  pgs_decline*transition

pam_pgs_strat <- mgcv::bam(
  formula = formula_pgs_strat,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_pgs_strat)

### plot survival curves
newData_pgs_strat <- make_newdata(ped_events, tend = sort(unique(tend)), transition=unique(transition), pgs_decline = c(0, 12, 24)) %>%
  group_by(transition, pgs_decline) %>%
  add_hazard(pam_pgs_strat) %>%
  add_cumu_hazard(pam_pgs_strat) %>%
  add_trans_prob(pam_pgs_strat) %>%
  ungroup() %>%
  mutate(transition = factor(transition, levels = c("0->1", "1->2", "2->3", "0->4", "1->4", "2->4")))

survCurves_pgs_strat <- ggplot(newData_pgs_strat, aes(x=tend, y=trans_prob)) +
  geom_line(aes(col=as.factor(pgs_decline))) +
  # geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper), alpha = 0.2) +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(min_age,80)) +
  ylim(c(0,1)) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(file.path(dir_out, "/figures/pam_pgs_strat_survCurves.png"), plot = survCurves_pgs_strat, width = 10, height = 10, units = "cm")

## genetics model (rs77924615) ----

### fit model
formula_genetics <- ped_status ~
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  sex*transition +
  rs77924615_A_G*transition +
  pgs_cross_594_umod*transition +
  diabetes*transition +
  smoking*transition +
  s(sc_BMI, by = transition) +
  s(sc_UACR, by = transition) +
  s(eGFRcrea, by = transition)

pam_genetics <- mgcv::bam(
  formula = formula_genetics,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_genetics)

png(file.path(dir_out, "/figures/pam_genetics_hazards.png"))
plot(pam_genetics, xlim=c(40,80), ylim=c(-5,5), page =1)
dev.off()

#### iteratively add risk factors and create table ----

# Calculate the sum of NAs for the specified columns
na_counts <- ped_events %>%
  summarise(across(
    c(tend, tend_onset, age_onset, tend_progression, age_progression, sex, transition, diabetes, sc_BMI, sc_UACR, eGFRcrea, pgs_cross_594_umod),
    ~ sum(is.na(.))
  ))

print(na_counts)

df_counts <- ped_events
# df_counts <- ped_events %>% filter(!is.na(pgs_cross_594_umod))
# df_counts <- ped_events %>% filter(!is.na(pgs_cross_594_umod), !is.na(sc_BMI), !is.na(sc_UACR))

idx_0 <- df_counts %>% filter(from == 0) %>% pull(id) %>% unique()
idx_1 <- df_counts %>% filter(from == 1) %>% pull(id) %>% unique()
idx_2 <- df_counts %>% filter(from == 2) %>% pull(id) %>% unique()

length(idx_0) # #at-risk for 0->1
df_counts %>% filter(transition == "0->1", ped_status == 1) %>% nrow() # #events for 0->1

length(idx_1) # #at-risk for 1->2
df_counts %>% filter(transition == "1->2", ped_status == 1) %>% nrow() # #events for 1->2

length(idx_2) # #at-risk for 2->3
df_counts %>% filter(transition == "2->3", ped_status == 1) %>% nrow() # #events for 2->3


# Helper to format p‐values in scientific notation: “mantissa x 10^exponent”
format_sci <- function(pval, digits = 2) {
  sci_str <- formatC(pval, format = "e", digits = digits)
  parts   <- strsplit(sci_str, "e")[[1]]
  mant    <- parts[1]
  expn    <- as.integer(parts[2])
  paste0(mant, " x 10^", expn)
}

# 1) Define the “always‐included” part of the formula:
base_terms <- "
  s(tend, by = transition_to_death) +
  s(tend_onset, by = transition_after_onset) +
  s(age_onset, by = transition_after_onset) +
  s(tend_progression, by = transition_after_progression) +
  s(age_progression, by = transition_after_progression) +
  sex * transition
"

# 2) For each iteration, list just the additional “rs77924615_A_G*transition”‐type terms:
extra_terms_list <- list(
  # 1st iteration: only “G × transition”
  "rs77924615_A_G * transition",

  # 2nd iteration: add PGS × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition"),

  # 3rd iteration: add Diabetes × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition"),

  # 4th iteration: add s(sc_BMI, by = transition) + s(sc_UACR, by = transition)
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, by = transition)",
    "s(sc_UACR, by = transition)"),

  # 5th iteration: add smoking × transition
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, by = transition)",
    "s(sc_UACR, by = transition)",
    "smoking * transition"),

  # 6th iteration: add s(eGFRcrea, by = transition)
  c("rs77924615_A_G * transition",
    "pgs_cross_594_umod * transition",
    "diabetes * transition",
    "s(sc_BMI, by = transition)",
    "s(sc_UACR, by = transition)",
    "smoking * transition",
    "s(eGFRcrea, by = transition)")

)

# 3) Friendly scenario names (one per iteration)
scenario_names <- c(
  "G only",
  "G + PGS",
  "G + PGS + Diabetes",
  "G + PGS + Diabetes + BMI + UACR",
  "G + PGS + Diabetes + BMI + UACR + Smoking",
  "G + PGS + Diabetes + BMI + UACR + Smoking + eGFR"
)

# 4) Pre‐allocate a list to store results
results <- vector("list", length(extra_terms_list))

# 5) Loop over each scenario, fit bam(), extract three coefficients, and store in a row
for (i in seq_along(extra_terms_list)) {
  extras <- extra_terms_list[[i]]

  # Build the right‐hand side: base_terms + extras
  rhs_terms <- paste0(base_terms, collapse = " ")
  rhs      <- paste( c(strsplit(base_terms, "\\+")[[1]] %>% trimws(), extras), collapse = " + " )

  # Create the full formula object
  fmla <- as.formula( paste("ped_status ~", rhs) )

  # Fit bam()
  pam_mod <- bam(
    formula  = fmla,
    data     = ped_events,
    family   = poisson(),
    offset   = offset,
    method   = "fREML",
    discrete = TRUE
  )

  # Summarize
  su <- summary(pam_mod)

  # Extract the three coefficients of interest (assumes exact naming in p.coeff)
  # NOTE: use the exact names: "rs77924615_A_G" (main),
  # "transition1->2:rs77924615_A_G", and "transition2->3:rs77924615_A_G"
  # If your factor‐coding yields a slightly different naming (e.g. "rs77924615_A_G:transition1->2"), adjust accordingly.
  coef_main  <- su$p.coeff["rs77924615_A_G"]
  se_main    <- su$se["rs77924615_A_G"]
  p_main     <- su$p.pv["rs77924615_A_G"]

  coef_t12   <- su$p.coeff["transition1->2:rs77924615_A_G"]
  se_t12     <- su$se["transition1->2:rs77924615_A_G"]
  p_t12      <- su$p.pv["transition1->2:rs77924615_A_G"]

  coef_t23   <- su$p.coeff["transition2->3:rs77924615_A_G"]
  se_t23     <- su$se["transition2->3:rs77924615_A_G"]
  p_t23      <- su$p.pv["transition2->3:rs77924615_A_G"]

  # Round the numeric estimates to two decimals, format p-values in scientific form
  coef_main_r  <- round(coef_main, 2)
  se_main_r    <- round(se_main,   2)
  p_main_fmt   <- format_sci(p_main, digits = 2)

  coef_t12_r   <- round(coef_t12, 2)
  se_t12_r     <- round(se_t12,   2)
  p_t12_fmt    <- format_sci(p_t12,  digits = 2)

  coef_t23_r   <- round(coef_t23, 2)
  se_t23_r     <- round(se_t23,   2)
  p_t23_fmt    <- format_sci(p_t23,  digits = 2)

  # Store in a one‐row tibble
  results[[i]] <- tibble(
    Scenario           = scenario_names[i],

    G_Coef             = coef_main_r,
    G_SE               = se_main_r,
    G_p_value          = p_main_fmt,

    `t1→2:G_Coef`      = coef_t12_r,
    `t1→2:G_SE`        = se_t12_r,
    `t1→2:G_p_value`   = p_t12_fmt,

    `t2→3:G_Coef`      = coef_t23_r,
    `t2→3:G_SE`        = se_t23_r,
    `t2→3:G_p_value`   = p_t23_fmt
  )
}

# 6) Bind all rows into the final 6×10 table
final_table <- bind_rows(results) %>%
  # 1) Convert all “_Coef” and “_SE” columns to character
  mutate(across(matches("_(Coef|SE)$"), ~ as.character(.))) %>%
  # 2) Now pivot to long form
  pivot_longer(
    cols          = -Scenario,
    names_to      = c("Effect", "Metric"),
    names_pattern = "(.+)_(Coef|SE|p_value)",
    values_to     = "Value"
  ) %>%
  # 3) Spread “Metric” back out to columns
  pivot_wider(
    names_from  = Metric,
    values_from = Value
  ) %>%
  # 4) Ensure the final column ordering
  select(Scenario, Effect, Coef, SE, p_value)

# 7) Print
print(final_table)

# 8) Save the final table to a CSV file
write.xlsx(final_table, file.path(dir_out, "genetics_model_results_originalCutoffs.xlsx"), rowNames = FALSE)

## genetics model (rs1047891_C_A) ----
formula_genetics_2 <- ped_status ~ s(tend) + s(tend_onset) + s(tend_progression) + transition + sex + age_onset + age_progression + rs1047891_C_A*transition

pam_genetics_2 <- mgcv::bam(
  formula = formula_genetics_2,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_genetics_2)

## null model (smooth onset age effect) ----

### fit model
formula_null_smoothAgeOnset <- ped_status ~ s(tend) + s(tend_onset) + s(tend_progression) + transition + sex + s(age_onset) + s(age_progression)

pam_null_smoothAgeOnset <- mgcv::bam(
  formula = formula_null_smoothAgeOnset,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_null_smoothAgeOnset)

png(file.path(dir_out, "/figures/pam_null_hazards_smoothAgeOnset.png"))
plot.gam(pam_null_smoothAgeOnset, select = 4)
dev.off()

png(file.path(dir_out, "/figures/pam_null_hazards_smoothAgeProgression.png"))
plot.gam(pam_null_smoothAgeOnset, select = 5)
dev.off()

## eGFR model ----

### fit model
formula_egfr <- ped_status ~ s(tend) + s(tend_onset) + s(tend_progression) + transition + sex + age_onset + age_progression + s(eGFRcrea, by = transition)

pam_egfr <- mgcv::bam(
  formula = formula_egfr,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_egfr)

png(file.path(dir_out, "/figures/pam_egfr_hazards_01.png"))
plot.gam(pam_egfr, select = 4, xlim = c(140, 60), ylim = c(-5,3))
dev.off()

png(file.path(dir_out, "/figures/pam_egfr_hazards_12.png"))
plot.gam(pam_egfr, select = 6, xlim = c(60, 30), ylim = c(-20,10))
dev.off()

png(file.path(dir_out, "/figures/pam_egfr_hazards_23.png"))
plot.gam(pam_egfr, select = 8, xlim=c(30,15), ylim = c(-5,3))
dev.off()

### plot cumulative hazards
png(file.path(dir_out, "/figures/pam_baseline_hazards.png"))
plot(pam_baseline, xlim = c(0,50), ylim = c(-20, 20), page=1)
dev.off()

### plot survival curves
ci <- FALSE
newData_baseline <- make_newdata(
    ped_events,
    tend = unique(tend),
    transition=unique(transition),
    # eGFRcrea = c(70, 90, 120)
    eGFRcrea = c(40,50),
    # age0 = quantile(age0, probs=c(0.10, 0.50, 0.90))

  ) %>%
  group_by(transition, eGFRcrea) %>%
  add_cumu_hazard(pam_baseline) %>%
  add_trans_prob(pam_baseline, ci = ci)

survCurves_baseline <- ggplot(newData_baseline, aes(x=tend, y=trans_prob)) +
  geom_line(aes(col=as.factor(eGFRcrea))) +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(30,80)) +
  ylim(c(0,1))
if (ci) {
  survCurves_baseline <- survCurves_baseline +
    geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper, fill = as.factor(age0)), alpha = 0.2)
}
ggsave(file.path(dir_out, "/figures/pam_baseline_survCurves.png"), plot = survCurves_baseline, width = 10, height = 10, units = "cm")

## genetics model ----

### fit model
ped_events$rs77924615_A_G_genotype <- ifelse(ped_events$rs77924615_A_G >= 1.5, 2, ifelse(ped_events$rs77924615_A_G >= 0.5, 1, 0))
formula_genetics <- ped_status ~ s(tend, by=interaction(transition,rs77924615_A_G_genotype)) + transition * rs77924615_A_G_genotype + s(eGFRcrea) + s(eGFRcrea, by=transition)

formula_genetics <- ped_status ~ s(tend, by=transition) + transition * rs77924615_A_G + s(eGFRcrea) + s(eGFRcrea, by=transition)

pam_genetics <- mgcv::bam(
  formula = formula_genetics,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_genetics)

### plot survival curves
ci <- FALSE
newData_genetics <- make_newdata(
    ped_events,
    # tend = unique(tend),
    tend = c(60),
    transition=unique(transition),
    rs77924615_A_G=c(0,1,2)
    # eGFRcrea = c(70, 90, 120)
    # eGFRcrea = c(40,50),
    # age0 = quantile(age0, probs=c(0.10, 0.50, 0.90))
  ) %>%
  group_by(transition, rs77924615_A_G) %>%
  add_cumu_hazard(pam_genetics) %>%
  add_trans_prob(pam_genetics, ci = ci)

survCurves_genetics <- ggplot(newData_genetics, aes(x=tend, y=trans_prob)) +
  geom_line(aes(col=as.factor(rs77924615_A_G))) +
  facet_wrap(~transition, ncol = 3) +
  xlim(c(30,80)) +
  ylim(c(0,1))
if (ci) {
  survCurves_genetics <- survCurves_genetics +
    geom_ribbon(aes(ymin = trans_lower, ymax = trans_upper, fill = as.factor(rs77924615_A_G)), alpha = 0.2)
}
ggsave(file.path(dir_out, "/figures/pam_genetics_survCurves.png"), plot = survCurves_genetics, width = 10, height = 10, units = "cm")

## PGS model ----

## tensor
# formula_pgs_decline <- ped_status ~ s(tend, by=transition) + te(pgs_decline, tend, by=transition) + s(eGFRcrea)
formula_pgs_decline <- ped_status ~ s(tend, by=transition) + pgs_decline*transition + s(eGFRcrea)

pam_pgs_decline <- mgcv::bam(
  formula = formula_pgs_decline,
  data = ped_events,
  family=poisson(),
  offset=offset,
  method="fREML",
  discrete = TRUE)

summary(pam_pgs_decline)

png(file.path(dir_out, "/figures/pam_pgs_decline.png"))
plot(pam_pgs_decline, xlim = c(35,80), ylim = c(-20, 20), page=1)
dev.off()

pam_extended_tensor <- gg_tensor(pam_pgs_decline, ci = TRUE) +
   xlab("age") +
   ylab("eGFR-at-baseline")
ggsave(file.path(dir_out, "/figures/pam_pgs_tensor_heatmap.png"), plot = pam_extended_tensor, width = 10, height = 10, units = "cm")

slice_age <- gg_slice(ped_events, pam_pgs_decline, term = "tend", pgs_decline = c(0,10,20), tend=35:80)
ggsave(file.path(dir_out, "/figures//pam_pgs_tensor_slice_age.png"), plot = slice_age, width = 10, height = 10, units = "cm")
slice_pgs <- gg_slice(ped_events, pam_pgs_decline, term = "pgs_decline", tend  = c(40,50,60,70), pgs_decline=0:24)
ggsave(file.path(dir_out, "/figures//pam_tensor_slice_eGFR.png"), plot = slice_pgs, width = 10, height = 10, units = "cm")

library(patchwork)

# Extract the unique transitions (assuming ped_events$transition is a factor or character)
transitions <- sort(unique(ped_events$transition))

# Create gg_slice() plots for the "tend" term (e.g., age effect) for each transition
slice_age_list <- lapply(transitions, function(tr) {
  gg_slice(as.data.frame(ped_events), pam_pgs_decline,
           term = "tend",
           # Fix pgs_decline at several values for comparison
           pgs_decline = c(0, 10, 20),
           # Evaluate predictions over a range of tend values
           tend = 35:80,
           # Fix transition to the current level so that the plot is stratified
           fixed = list(transition = tr)) +
    ggtitle(paste("Transition:", tr))
})

# Combine the individual tend plots into one figure (one per row; adjust ncol/nrow as desired)
combined_slice_age <- wrap_plots(slice_age_list, ncol = 1)
combined_slice_age

# Now, similarly, create gg_slice() plots for the "pgs_decline" term for each transition
slice_pgs_list <- lapply(transitions, function(tr) {
  gg_slice(ped_events, pam_pgs_decline,
           term = "pgs_decline",
           # Evaluate predictions for various tend values
           tend = c(40, 50, 60, 70),
           # Sweep pgs_decline over a range of values
           pgs_decline = 0:24,
           fixed = list(transition = tr)) +
    ggtitle(paste("Transition:", tr))
})

# Combine these plots into one figure (again, adjust layout as needed)
combined_slice_pgs <- wrap_plots(slice_pgs_list, ncol = 1)
combined_slice_pgs
