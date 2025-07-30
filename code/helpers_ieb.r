library(data.table)
library(dplyr)
library(tidyr)
library(rlang) # for sim_pexp
library(Formula) # for sim_pexp
library(lazyeval) # for sim_pexp
library(purrr) # for sim_pexp
library(survival)
library(mgcv)
library(pammtools)
library(flexsurv)
library(ggplot2)
library(scales)
library(patchwork)

source("/nvmetmp/wis37138/msm_kidneyFunction/code/helpers_sim.r")

# data prep ----

wrapper_sim <- function(
  job,
  data,
  dist_x1 = "rbinom(n, 1, 0.5)",
  dist_x2 = "rbinom(n, 1, 0.5)",
  beta_vals,
  formulas_dgp,
  terminal_states,
  cut,
  n = 2500,
  round = 2,
  cens_type = c("none", "right", "interval", "left"),
  cens_dist = c("weibull", "exponential", "lognormal", "uniform"),
  cens_params = NULL
  ){

  df <- data.frame(
    id          = seq_len(n),
    from        = 0L,
    t           = 0
  )

  vars <- get_x_vars(formulas_dgp)
  dist_list <- list(x1 = dist_x1, x2 = dist_x2)
  for (var in vars) {
    if (!var %in% names(dist_list)) {
      stop("No distribution provided for variable: ", var)
    }
    df[[var]] <- eval(parse(text = dist_list[[var]]))
  }

  beta_vals_unlisted <- unlist(beta_vals)

  if(dist_x1 == "rnorm(n, 0, 5)" &
     dist_x2 == "rnorm(n, 0, 5)" &
     beta_vals_unlisted["beta_1_01"] == 0.5) {
      beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
      beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 6.0
     }

  formulas_dgp <- lapply(formulas_dgp, function(fli) {
    old_f  <- fli$formula          # e.g. ~ f_0(tend) + beta_1_01*x1 + …
    rhs_old <- old_f[[2]]          # extract the RHS expression
    rhs_new <- bake_expr(rhs_old, beta_vals_unlisted)
    # rebuild a formula "~ <rhs_new>" in the global env so f_0/g_0 etc still resolve
    fli$formula <- as.formula(call("~", rhs_new), env = globalenv())
    fli
  })

  events <- sim_pexp_msm(formulas_dgp, df, cut, terminal_states, round = round, add_counterfactuals = FALSE)

  if(cens_type != "none") {
    events <- events %>%
      add_censoring(type = cens_type, distribution = cens_dist, parameters = cens_params, round = round)
  }

  ped <- as_ped_multistate(
    data       = events %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

  ped <- ped %>%
    add_timeScales() %>%
    add_transVars()

  out <- list(
    events = events,
    ped = ped,
    formulas_dgp  = formulas_dgp,
    beta_vals = beta_vals_unlisted
  )
saveRDS(out, "out.rds")
  return(out)
}


wrapper_msm <- function(
  data,
  job,
  instance,
  formula) {

  ped <- instance$ped

  ped$trans_after_1 <- factor(ped$trans_after_1,
                              levels = c("none", "1->2", "1->3"),
                              ordered = TRUE)

  mod <- bam(as.formula(formula)
                      , data = ped
                      , family=poisson()
                      , offset=offset
                      , discrete = T
                      , method = "fREML"
  )
  summary_mod <- summary(mod)
  vcov_mod <- vcov(mod)

  formulas_dgp <- instance$formulas_dgp

  # beta_vars <- formulas_dgp %>%
  #   map(~ {
  #     vars <- all.vars(.x$formula)
  #     vars[grepl("^beta_[12]_(01|12)$", vars)]
  #   }) %>%
  #   unlist() %>%
  #   unique() %>%
  #   sort()

  # beta_true <- beta_vars %>%
  #   map(~ eval(parse(text = .x), envir = ped)) %>%
  #   set_names(beta_vars)

  beta_true <- as.list(instance$beta_vals)

  # Create a dataframe with coefficient information
  coef_df <- data.frame(
    transition = c("onset", "progression_int", "progression"),
    coefficient = c(summary_mod$p.coeff["x1"],
                   summary_mod$p.coeff["transition1->2:x1"],
                   summary_mod$p.coeff["x1"] + summary_mod$p.coeff["transition1->2:x1"]),
    std_error = c(summary_mod$se["x1"],
                 summary_mod$se["transition1->2:x1"],
                 sqrt(vcov_mod["x1", "x1"] + vcov_mod["transition1->2:x1", "transition1->2:x1"] +
                        2 * vcov_mod["x1", "transition1->2:x1"])),
    p_value = c(summary_mod$p.pv["x1"],
               summary_mod$p.pv["transition1->2:x1"],
               NA),
    stringsAsFactors = FALSE
  )

  coef_df[coef_df$transition == "progression", "p_value"] <-
    2 * (1 - pnorm(abs(coef_df[coef_df$transition == "progression", "coefficient"] /
                      coef_df[coef_df$transition == "progression", "std_error"])))

  # Round numeric values for better readability
  coef_df <- coef_df %>%
    mutate(
      coefficient = round(coefficient, 4),
      std_error = round(std_error, 4),
      p_value = round(p_value, 4),
      beta_1 = c(beta_true$beta_1_01, beta_true$beta_1_12 - beta_true$beta_1_01, beta_true$beta_1_12),
      beta_2 = c(beta_true$beta_2_01, beta_true$beta_2_12 - beta_true$beta_2_01, beta_true$beta_2_12),
      coverage = ifelse(
        coefficient - 1.96 * std_error <= beta_1 &
        coefficient + 1.96 * std_error >= beta_1,
        1, 0
      ),
      bias = coefficient - beta_1
    )

  return(coef_df)
}


wrapper_cor <- function(
  data,
  job,
  instance) {

  df <- instance$events

  vars <- get_x_vars(instance$formulas_dgp)
  if (length(vars) == 0) {
    stop("No covariates found in the formulas.")
  }
  if (!all(vars %in% colnames(df))) {
    missing_vars <- vars[!vars %in% colnames(df)]
    stop("The following variables are not found in the data: ",
         paste(missing_vars, collapse = ", "))
  }

  cor_df <- calc_statewise_correlation(
    data = df,
    from_states = unique(df$from),
    vars = vars
  )

  return(cor_df)
}


# Recursively walk an R expression and replace any symbol in `vals`
# with its numeric value.
bake_expr <- function(expr, vals) {
  if (is.name(expr)) {
    nm <- as.character(expr)
    if (nm %in% names(vals)) {
      return(vals[[nm]])        # replace symbol with number
    } else {
      return(expr)              # keep other symbols as-is
    }
  }
  if (!is.call(expr)) {
    return(expr)                # e.g. constants, literals
  }
  # it's a call: recurse over each element
  elts <- as.list(expr)
  elts <- lapply(elts, bake_expr, vals)
  return(as.call(elts))
}


calc_statewise_correlation <- function(data, from_states, vars) {
  # Check if vars contains valid column names
  if (!all(vars %in% colnames(data))) {
    missing_vars <- vars[!vars %in% colnames(data)]
    stop("The following variables are not found in the data: ",
         paste(missing_vars, collapse = ", "))
  }

  # Initialize an empty dataframe to store results
  results_df <- data.frame(
    from = character(),
    var1 = character(),
    var2 = character(),
    rho = numeric(),
    p = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through each state in from_states
  for (state in from_states) {
    # Filter data for current state
    state_data <- data[data$from == state, , drop = FALSE]

    # If no data for this state, skip
    if (nrow(state_data) == 0) {
      next
    }

    # Loop through all pairs of variables
    for (i in 1:(length(vars) - 1)) {
      for (j in (i + 1):length(vars)) {
        var1 <- vars[i]
        var2 <- vars[j]

        # Calculate correlation and p-value using cor.test
        test_result <- cor.test(state_data[[var1]], state_data[[var2]],
                                use = "pairwise.complete.obs")

        # Add row to results dataframe
        results_df <- rbind(results_df, data.frame(
          from = as.character(state),
          var1 = var1,
          var2 = var2,
          rho = test_result$estimate,
          p = test_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  return(results_df)
}


plot_cor <- function(df, var_name = "rho") {
  df <- df %>%
    mutate(
      from    = factor(from, levels = c("0", "1")),
      problem = recode(problem, !!!problem_labs)
    )

  true_val <- 0

  ggplot(df,
         aes(
           x     = from,
           y     = rho,
           fill  = from,                        # <— map fill to state
           group = interaction(from, algorithm) # keeps separate boxes per state×algo
         )) +

    geom_hline(
      yintercept = true_val,
      colour     = "orange",
      linewidth  = 1.3,
      linetype   = "solid"
    ) +

    geom_boxplot(
      outlier.shape = NA,
      position      = position_dodge(width = 0.75)
    ) +

    # manually set state‐colours: 0 = green, 1 = red
    scale_fill_manual(
      name   = "State",
      values = c("0" = "green", "1" = "red"),
      labels = c("0", "1")
    ) +

    facet_wrap(~ problem, nrow = 1) +

    scale_x_discrete(
      name   = "State",
      labels = c("0", "1")
    ) +

    # Added scale_y_continuous for requested limits/breaks
    scale_y_continuous(
      limits = c(-0.51, 0.05),
      breaks = seq(-0.50, 0.0, by = 0.10),
      labels = number_format(accuracy = 0.01)
    ) +

    labs(
      y     = var_name
    ) +

    theme_bw(base_size = 14) +
    theme(
      axis.text.x     = element_text(angle = 0),
      panel.spacing.x = unit(1.2, "lines")
    )
}


plot_coefs <- function(df, var_name = "coefficient") {
  df2 <- df %>%
    mutate(
      algo       = case_when(
        grepl("^stratified", formula_name) ~ "stratified",
        grepl("^timeScales", formula_name) ~ "timeScales",
        TRUE ~ formula_name
      ),
      ieb        = if_else(grepl("_ieb$", formula_name), "IEB", "Base"),
      problem    = recode(problem, !!!problem_labs),
      transition = factor(transition, levels = unique(transition))
    )

  df_hlines <- df2 %>%
    distinct(transition, problem, beta_1)

  # Split data by transition
  p_list <- lapply(levels(df2$transition), function(tran) {
    df_sub <- df2 %>% filter(transition == tran)
    df_hlines_sub <- df_hlines %>% filter(transition == tran)
    # Set y-limits and breaks based on transition
    if (tran %in% c("onset", "progression")) {
      ylims <- c(0.05, 0.60)
      ybreaks <- seq(0.1, 0.6, by = 0.1)
    } else if (tran == "progression_int") {
      ylims <- c(-0.20, 0.08)
      ybreaks <- seq(-0.2, 0.1, by = 0.1)
    } else {
      ylims <- NULL; ybreaks <- waiver()
    }

    ggplot(df_sub, aes(
      x     = algo,
      y     = coefficient,
      fill  = ieb,
      color = ieb,
      group = interaction(algo, ieb)
    )) +
      geom_hline(
        data        = df_hlines_sub,
        aes(yintercept = beta_1),
        colour      = "orange",
        linewidth   = 1.3,
        linetype    = "solid"
      ) +
      geom_boxplot(
        outlier.shape = NA,
        position      = position_dodge(width = 0.75)
      ) +
      scale_fill_brewer(
        palette = "Dark2",
        name    = "Model",
        labels  = c("Base", "IEB")
      ) +
      scale_color_brewer(
        palette = "Dark2",
        guide   = "none"
      ) +
      facet_grid(cols = vars(problem)) +
      scale_y_continuous(
        limits = ylims,
        breaks = ybreaks,
        labels = number_format(accuracy = 0.1)
      ) +
      labs(
        x = "Algorithm",
        y = var_name,
        title = tran
      ) +
      theme_bw(base_size = 14) +
      theme(
        axis.text.x   = element_text(angle = 0),
        panel.spacing = unit(1.2, "lines")
      )
  })

  # Combine plots vertically
  p <- p_list[[1]] / p_list[[2]] / p_list[[3]]

  return(p)
}


plot_coverage <- function(df) {
  # 1) derive “algo” and flag IEB, recode problem, order transitions
  df2 <- df %>%
    mutate(
      algo       = case_when(
                     grepl("^stratified", formula_name) ~ "stratified",
                     grepl("^timeScales", formula_name) ~ "timeScales",
                     TRUE ~ formula_name
                   ),
      ieb        = if_else(grepl("_ieb$", formula_name), "IEB", "Base"),
      problem    = recode(problem, !!!problem_labs),
      transition = factor(transition, levels = unique(transition))
    )

  # 2) summarise & run binom.test per group
  df_sum <- df2 %>%
    group_by(transition, problem, algo, ieb) %>%
    summarise(
      mean_cov = mean(coverage, na.rm = TRUE),
      x        = sum(coverage, na.rm = TRUE),
      n        = n(),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      bt        = list(binom.test(x, n, p = 0.5)),
      ci_lower  = bt$conf.int[1],
      ci_upper  = bt$conf.int[2]
    ) %>%
    ungroup()

  # 3) plot
  ggplot(df_sum, aes(x = algo, y = mean_cov, fill = ieb)) +
    # bars
    geom_col(position = position_dodge(width = 0.75)) +
    # black CI bars
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                  position = position_dodge(width = 0.75),
                  width    = 0.2,
                  colour   = "black") +
    # dashed red 95% target line
    geom_hline(yintercept = 0.95,
               colour      = "red",
               linetype    = "dashed",
               linewidth   = 1) +
    # percent scale with 95% labelled
    scale_y_continuous(
      name   = "Coverage",
      breaks = c(seq(0, 1, by = 0.2), 0.95),
      labels = percent_format(accuracy = 1)
    ) +
    # fill palette
    scale_fill_brewer(
      palette = "Dark2",
      name    = "Model",
      labels  = c("Base", "IEB")
    ) +
    # facets
    facet_grid(rows = vars(transition),
               cols = vars(problem)) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1),
      panel.spacing = unit(1.2, "lines")
    )
}
