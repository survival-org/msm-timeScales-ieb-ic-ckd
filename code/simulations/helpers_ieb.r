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

  if(
    dist_x2 == "rnorm(n, 0, 5)" &
    beta_vals_unlisted["beta_1_01"] == 0.5) {

    if(dist_x1 == "rnorm(n, 0, 5)") {
      beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
      beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 6.0
     } else if(dist_x1 == "rnorm(n, 0, 1)") {
      beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
      beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 3.0
    } else if(dist_x1 == "rbinom(n, 1, 0.5)") {
      beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
      beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 3.0
    }
  } else if(
    dist_x1 == "rnorm(n, 0, 5)" &
    dist_x2 != "rnorm(n, 0, 5)" &
    beta_vals_unlisted["beta_1_01"] == 0.5) {

      if(dist_x2 == "rnorm(n, 0, 1)") {
        beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
        beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 3.0
      } else if(dist_x2 == "rbinom(n, 1, 0.5)") {
        beta_vals_unlisted["beta_0_01"] <- beta_vals_unlisted["beta_0_01"] - 2.0
        beta_vals_unlisted["beta_0_12"] <- beta_vals_unlisted["beta_0_12"] - 3.0
      }

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


summarize_bias_results <- function(data) {
  # Group by specified variables, now including beta_scenario
  summary_data <- data %>%
    dplyr::group_by(transition, problem, formula_name, dist_x1_x2, beta_scenario) %>%
    dplyr::summarise(
      coef_mean = mean(coefficient, na.rm = TRUE),
      bias_mean = mean(bias, na.rm = TRUE),
      # Calculate empirical 95% CI for the bias
      bias_lower = quantile(bias, 0.025, na.rm = TRUE),
      bias_upper = quantile(bias, 0.975, na.rm = TRUE),
      coverage_mean = mean(coverage, na.rm = TRUE),
      .groups = "drop" # Drop grouping for subsequent steps
    )

  # Process columns and filter rows
  processed_summary <- summary_data %>%
    # Filter for only the specified transitions
    dplyr::filter(transition %in% c("onset", "progression")) %>%
    # Split dist_x1_x2 into two columns
    tidyr::separate(dist_x1_x2, into = c("dist_x1", "dist_x2"), sep = "_", remove = FALSE) %>%
    # Recode distribution names into a readable format
    dplyr::mutate(
      dist_x1 = dplyr::case_when(
        dist_x1 == "bernoulli0.5" ~ "Bernoulli(0.5)",
        dist_x1 == "normal1" ~ "N(0,1)",
        dist_x1 == "normal5" ~ "N(0,5)",
        TRUE ~ as.character(dist_x1)
      ),
      dist_x2 = dplyr::case_when(
        dist_x2 == "bernoulli0.5" ~ "Bernoulli(0.5)",
        dist_x2 == "normal1" ~ "N(0,1)",
        dist_x2 == "normal5" ~ "N(0,5)",
        TRUE ~ as.character(dist_x2)
      )
    )

  return(processed_summary)
}


plot_cor <- function(df, var_name = "rho", font_size = 14, y_limits = NULL, y_breaks = NULL) {
  df <- df %>%
    mutate(
      from    = factor(from, levels = c("0", "1")),
      problem = recode(problem, !!!problem_labs),
      problem = factor(problem, levels = c("SSTS DGP", "MTS DGP"))
    )

  true_val <- 0

  ggplot(df,
         aes(
           x     = from,
           y     = rho,
           fill  = from,
           color = from,
           group = interaction(from, algorithm)
         )) +

    geom_hline(
      yintercept = true_val,
      colour     = "orange",
      linewidth  = 1.3,
      linetype   = "solid"
    ) +

    geom_boxplot(
      # outlier.shape = NA,
      position      = position_dodge(width = 0.75)
    ) +

    geom_boxplot(
      position      = position_dodge(width = 0.75),
      color         = "black",
      fill          = NA,
      coef          = 0,
      outlier.shape = NA
    ) +

    scale_fill_manual(
      name   = "State",
      values = c("0" = "#2c65cf", "1" = "darkgreen"),
      labels = c("0", "1")
    ) +

    scale_color_manual(
      name   = "State",
      values = c("0" = "#2c65cf", "1" = "darkgreen"),
      guide  = "none"
    ) +

    facet_wrap(~ problem, nrow = 1) +

    scale_x_discrete(
      labels = c("0", "1")
    ) +

    scale_y_continuous(
      limits = y_limits,
      breaks = y_breaks,
      labels = number_format(accuracy = 0.01)
    ) +

    labs(
      y     = "Correlation",
      x     = "State",
      title = NULL
    ) +

    theme_bw(base_size = font_size) +
    theme(
      axis.text.x     = element_text(size = font_size, angle = 0),
      axis.text.y     = element_text(size = font_size),
      axis.title.x    = element_text(size = font_size),
      axis.title.y    = element_text(size = font_size),
      strip.text      = element_text(size = font_size),
      plot.title      = element_text(size = font_size, hjust = 0.5),
      legend.text     = element_text(size = font_size),
      legend.title    = element_text(size = font_size),
      panel.spacing.x = unit(1.2, "lines"),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

plot_coefs <- function(df, var_name = "coefficient", drop_int = TRUE, font_size = 14, y_limits = NULL, y_breaks = NULL) {

  if(drop_int) {
    df <- df %>%
      filter(transition != "progression_int")
  }

  df2 <- df %>%
    mutate(
      algo       = case_when(
        grepl("^stratified", formula_name) ~ "SSTS PAM",
        grepl("^timeScales", formula_name) ~ "MTS PAM",
        TRUE ~ formula_name
      ),
      algo       = factor(algo, levels = c("SSTS PAM", "MTS PAM")),
      # ieb        = if_else(grepl("_ieb$", formula_name), "IEB", "Oracle"),
      # ieb        = factor(ieb, levels = c("Oracle", "IEB")),
      ieb        = if_else(grepl("_ieb$", formula_name), "$x_1$-only Model", "Full Model"),
      ieb        = factor(ieb, levels = c("$x_1$-only Model", "Full Model")),
      problem    = recode(problem, !!!problem_labs),
      problem = factor(problem, levels = c("SSTS DGP", "MTS DGP")),
      transition = recode(transition,
              onset = "Onset",
              progression_int = "Progression Diff.",
              progression = "Progression"),
      transition = factor(transition, levels = c("Onset", "Progression Diff.", "Progression"))
    )

  df_hlines <- df2 %>%
    distinct(transition, problem, beta_1)

  p_list <- lapply(sort(unique(df2$transition)), function(tran) {
    tran <- as.character(tran)
    df_sub <- df2 %>% filter(transition == tran)
    df_hlines_sub <- df_hlines %>% filter(transition == tran)

    if (tran %in% c("Onset", "Progression")) {
      ylims <- y_limits[[1]]
      ybreaks <- y_breaks[[1]]

      if(tran == "Onset") {
        xlab <- ",0→1"
      } else if(tran == "Progression") {
        xlab <- ",1→2"
      }

    } else if (tran == "Progression Diff.") {
      ylims <- y_limits[[2]]
      ybreaks <- y_breaks[[2]]
      xlab <- ",1→2 (int)"
    } else {
      ylims <- NULL; ybreaks <- waiver(); xlab <- NULL
    }

    ggplot(df_sub, aes(
      x     = algo, y     = coefficient, fill  = ieb,
      color = ieb, group = interaction(algo, ieb)
    )) +
      geom_hline(
        data        = df_hlines_sub, aes(yintercept = beta_1),
        colour      = "orange", linewidth   = 1.3, linetype    = "solid"
      ) +
      geom_boxplot(
        # outlier.shape = NA,
        position = position_dodge(width = 0.75)) +
      geom_boxplot(
        position      = position_dodge(width = 0.75),
        color         = "black",
        fill          = NA,
        coef          = 0,
        outlier.shape = NA
      ) +
      scale_fill_manual(name = NULL, values = c("Full Model" = "#009E73", "$x_1$-only Model" = "#999999"), labels = c(expression(paste(x[1], "-only Model")), "Full Model")) +
      scale_color_manual(name = NULL, values = c("Full Model" = "#009E73", "$x_1$-only Model" = "#999999"), guide = "none") +
      facet_grid(cols = vars(problem)) +
      scale_y_continuous(limits = ylims, breaks = ybreaks, labels = number_format(accuracy = 0.1)) +
      # labs(y = as.expression(bquote(.(tran) ~ "(" * hat(beta)[x[.(paste0("1", xlab))]] * ")"))) +
      labs(y = tran) +
      theme_bw(base_size = font_size) +
      theme(
        axis.text.x   = element_text(size = font_size, angle = 0),
        axis.text.y   = element_text(size = font_size),
        axis.title.x  = element_blank(),
        axis.title.y  = element_text(size = font_size),
        strip.text    = element_text(size = font_size),
        plot.title    = element_text(size = font_size, hjust = 0.5),
        legend.text   = element_text(size = font_size),
        legend.title  = element_text(size = font_size),
        panel.spacing = unit(1.2, "lines"),
        legend.position = "bottom", legend.direction = "horizontal"
      )
  })

  return(p_list)
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
      ieb        = if_else(grepl("_ieb$", formula_name), "IEB", "Oracle"),
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
    # points instead of bars (filled circles with black outline)
    geom_point(
      position = position_dodge(width = 0.75),
      shape    = 21,
      size     = 3,
      colour   = "black"
    ) +
    # black CI bars (same dodge as points)
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      position = position_dodge(width = 0.75),
      width    = 0.2,
      colour   = "black"
    ) +
    # dashed red 95% target line
    geom_hline(yintercept = 0.95, colour = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(
      name   = "Coverage",
      breaks = c(seq(0, 1, by = 0.2), 0.95),
      labels = scales::percent_format(accuracy = 1)
    ) +
    scale_fill_brewer(palette = "Dark2", name = "Model", labels = c("Oracle", "IEB")) +
    facet_grid(rows = vars(transition), cols = vars(problem)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(1.2, "lines"))

}

process_correlation_summary <- function(data) {
  data %>%
    separate(dist_x1_x2, into = c("dist_x1", "dist_x2"), sep = "_", remove = FALSE) %>%
    mutate(
      dist_x1 = case_when(
        dist_x1 == "bernoulli0.5" ~ "Ber(0.5)",
        dist_x1 == "normal1" ~ "N(0,1)",
        dist_x1 == "normal5" ~ "N(0,5)",
        TRUE ~ as.character(dist_x1)
      ),
      dist_x2 = case_when(
        dist_x2 == "bernoulli0.5" ~ "Ber(0.5)",
        dist_x2 == "normal1" ~ "N(0,1)",
        dist_x2 == "normal5" ~ "N(0,5)",
        TRUE ~ as.character(dist_x2)
      )
    )
}


convert_to_latex_effect_sizes <- function(betas, size_names = NULL) {
  latex_table <- "\\begin{table}[htbp]\n"
  latex_table <- paste0(latex_table, "\\centering\n")
  latex_table <- paste0(latex_table, "\\caption{\\captioniebeffectsizes}\n")
  latex_table <- paste0(latex_table, "\\label{tab:sim-ieb-effect-sizes}\n")

  col_defs <- "{p{0.38\\textwidth} >{\\centering\\arraybackslash}p{0.14\\textwidth} >{\\centering\\arraybackslash}p{0.14\\textwidth} >{\\centering\\arraybackslash}p{0.14\\textwidth}}"
  latex_table <- paste0(latex_table, "\\begin{tabular}", col_defs, "\n")
  latex_table <- paste0(latex_table, "\\toprule\n")

  latex_table <- paste0(latex_table, "\\multicolumn{1}{c}{} & \\multicolumn{3}{c}{Effect Sizes} \\\\\n")
  latex_table <- paste0(latex_table, "\\cmidrule(lr){2-4}\n")

  if (is.null(size_names)) {
    capitalize <- function(s) {
      paste0(toupper(substring(s, 1, 1)), substring(s, 2))
    }
    final_headers <- sapply(names(betas), capitalize)
  } else {
    final_headers <- size_names
  }

  header_cols <- paste0("\\multicolumn{1}{c}{", final_headers, "}", collapse = " & ")
  header_line <- paste0("\\multicolumn{1}{l}{Coefficient} & ", header_cols, " \\\\\n")
  latex_table <- paste0(latex_table, header_line)

  latex_table <- paste0(latex_table, "\\midrule\n")

  format_beta_name <- function(name) {
    name <- gsub("beta_0_01", "$\\\\beta_{x_0,0\\\\rightarrow 1}$", name)
    name <- gsub("beta_0_03", "$\\\\beta_{x_0,0\\\\rightarrow 3}$", name)
    name <- gsub("beta_0_12", "$\\\\beta_{x_0,1\\\\rightarrow 2}$", name)
    name <- gsub("beta_0_13", "$\\\\beta_{x_0,1\\\\rightarrow 3}$", name)
    name <- gsub("beta_1_01", "$\\\\beta_{x_1,0\\\\rightarrow 1}$", name)
    name <- gsub("beta_1_03", "$\\\\beta_{x_1,0\\\\rightarrow 3}$", name)
    name <- gsub("beta_1_12", "$\\\\beta_{x_1,1\\\\rightarrow 2}$", name)
    name <- gsub("beta_1_13", "$\\\\beta_{x_1,1\\\\rightarrow 3}$", name)
    name <- gsub("beta_2_01", "$\\\\beta_{x_2,0\\\\rightarrow 1}$", name)
    name <- gsub("beta_2_03", "$\\\\beta_{x_2,0\\\\rightarrow 3}$", name)
    name <- gsub("beta_2_12", "$\\\\beta_{x_2,1\\\\rightarrow 2}$", name)
    name <- gsub("beta_2_13", "$\\\\beta_{x_2,1\\\\rightarrow 3}$", name)
    return(name)
  }

  coefficient_names <- names(betas[[1]])

  latex_table <- paste0(latex_table, "\\multicolumn{4}{l}{\\textbf{Transition-specific intercepts}} \\\\\n")
  for (i in 1:4) {
    beta_name <- coefficient_names[i]
    formatted_name <- format_beta_name(beta_name)

    val1 <- sprintf("%.1f", betas[[1]][[beta_name]])
    val2 <- sprintf("%.1f", betas[[2]][[beta_name]])
    val3 <- sprintf("%.1f", betas[[3]][[beta_name]])

    latex_table <- paste0(latex_table, formatted_name, " & ", val1, " & ", val2, " & ", val3, " \\\\\n")
  }

  latex_table <- paste0(latex_table, "\\midrule\n")

  latex_table <- paste0(latex_table, "\\multicolumn{4}{l}{\\textbf{Transition-specific effect sizes of risk factor }$x_1$} \\\\\n")
  for (i in 5:8) {
    beta_name <- coefficient_names[i]
    formatted_name <- format_beta_name(beta_name)

    val1 <- sprintf("%.1f", betas[[1]][[beta_name]])
    val2 <- sprintf("%.1f", betas[[2]][[beta_name]])
    val3 <- sprintf("%.1f", betas[[3]][[beta_name]])

    latex_table <- paste0(latex_table, formatted_name, " & ", val1, " & ", val2, " & ", val3, " \\\\\n")
  }

  latex_table <- paste0(latex_table, "\\midrule\n")

  latex_table <- paste0(latex_table, "\\multicolumn{4}{l}{\\textbf{Transition-specific effect sizes of risk factor }$x_2$} \\\\\n")
  for (i in 9:12) {
    beta_name <- coefficient_names[i]
    formatted_name <- format_beta_name(beta_name)

    val1 <- sprintf("%.1f", betas[[1]][[beta_name]])
    val2 <- sprintf("%.1f", betas[[2]][[beta_name]])
    val3 <- sprintf("%.1f", betas[[3]][[beta_name]])

    latex_table <- paste0(latex_table, formatted_name, " & ", val1, " & ", val2, " & ", val3, " \\\\\n")
  }

  latex_table <- paste0(latex_table, "\\bottomrule\n")
  latex_table <- paste0(latex_table, "\\end{tabular}\n")
  latex_table <- paste0(latex_table, "\\end{table}\n")

  return(latex_table)
}


convert_to_latex_cor <- function(data, round_digits = 2, caption = "\\captioniebcor", label = "tab:sim-ieb-cor", colsep = "2pt") {

  wide_data <- data %>%
    dplyr::mutate(
      dist_x1 = factor(dist_x1, levels = c("Ber(0.5)", "N(0,1)", "N(0,5)")),
      dist_x2 = factor(dist_x2, levels = c("Ber(0.5)", "N(0,1)", "N(0,5)"))
    ) %>%
    dplyr::arrange(dist_x1, dist_x2, beta_scenario) %>%
    dplyr::mutate(dplyr::across(c(rho_mean, rho_lower, rho_upper), ~ round(.x, round_digits))) %>%
    dplyr::mutate(cell_content = sprintf("\\makecell{%.*f \\\\ (%.*f; %.*f)}",
                                         round_digits, rho_mean,
                                         round_digits, rho_lower,
                                         round_digits, rho_upper)) %>%
    tidyr::unite(col_name, problem, from) %>%
    dplyr::select(dist_x1, dist_x2, beta_scenario, col_name, cell_content) %>%
    tidyr::pivot_wider(names_from = col_name, values_from = cell_content, values_fill = "-")

  col_order <- c("dist_x1", "dist_x2", "beta_scenario",
                 "sim_stratified_0", "sim_stratified_1",
                 "sim_timeScales_0", "sim_timeScales_1")

  missing_cols <- setdiff(col_order, names(wide_data))
  if (length(missing_cols) > 0) {
      for (col in missing_cols) {
        wide_data[[col]] <- "-"
      }
  }
  wide_data <- wide_data[, col_order]

  format_dist_name <- function(name) {
      name <- stringr::str_replace_all(name, c("Bernoulli" = "\\\\text{Bernoulli}", "N" = "N", "Ber" = "\\\\text{Ber}"))
      sprintf("\\(%s\\)", name)
  }

  body_rows <- apply(wide_data, 1, function(row) {
    dist1 <- format_dist_name(row["dist_x1"])
    dist2 <- format_dist_name(row["dist_x2"])
    beta_s <- row["beta_scenario"]
    paste(c(dist1, dist2, beta_s, row[4:7]), collapse = " & ")
  })

  body <- paste(body_rows, collapse = " \\\\\n")

  # Assemble the final LaTeX string with explicit line breaks and corrected header
  final_latex_code <- paste(
    "\\begin{table}[ht!]\n",
    "\\centering\n",
    paste0("\\setlength{\\tabcolsep}{", colsep, "}\n"),
    paste0("\\caption{", caption, "}\n"),
    paste0("\\label{", label, "}\n"),
    "\\footnotesize\n",
    "\\begin{tabular}{llclcccc}\n",
    "\\toprule\n",
    " & & \\textbf{Effect} & \\multicolumn{2}{c}{\\textbf{STSS DGP}} & \\multicolumn{2}{c}{\\textbf{MTS DGP}} \\\\\n", # "Effect" on its own line
    "\\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n",
    "\\textbf{Dist. \\(x_1\\)} & \\textbf{Dist. \\(x_2\\)} & \\textbf{Sizes} & \\textbf{State 0} & \\textbf{State 1} & \\textbf{State 0} & \\textbf{State 1} \\\\\n", # "Sizes" aligned with Dist. x1 and x2
    "\\midrule\n",
    body,
    "\\\\\n",
    "\\bottomrule\n",
    "\\normalsize\n",
    "\\end{tabular}\n",
    "\\end{table}",
    collapse = ""
  )

  return(invisible(final_latex_code))
}


create_csv_summary <- function(data) {
  # This function requires the dplyr, tidyr, and stringr packages.

  # 1. Summarize the raw simulation data
  summary_data <- data %>%
    dplyr::group_by(transition, problem, formula_name, dist_x1_x2, beta_scenario) %>%
    dplyr::summarise(
      # NEW: Keep the true coefficient value for the group
      coef_true = dplyr::first(beta_1),

      # Bias Summaries
      bias_mean = mean(bias, na.rm = TRUE),
      bias_lower = quantile(bias, 0.025, na.rm = TRUE),
      bias_upper = quantile(bias, 0.975, na.rm = TRUE),

      # Coefficient Summaries
      coef_mean = mean(coefficient, na.rm = TRUE),
      coef_lower = quantile(coefficient, 0.025, na.rm = TRUE),
      coef_upper = quantile(coefficient, 0.975, na.rm = TRUE),

      # Coverage Summaries (with binomial confidence interval)
      coverage_mean = mean(coverage, na.rm = TRUE),
      coverage_ci = list(stats::binom.test(sum(coverage), n())$conf.int),

      .groups = "drop"
    ) %>%
    tidyr::unnest_wider(coverage_ci, names_sep = "_") %>%
    dplyr::rename(coverage_lower = coverage_ci_1, coverage_upper = coverage_ci_2)

  # 2. Clean, transform, sort, and rename columns for presentation
  csv_output <- summary_data %>%
    # Create descriptive factor columns
    tidyr::separate(dist_x1_x2, into = c("dist_x1_name", "dist_x2_name"), sep = "_", remove = FALSE) %>%
    dplyr::mutate(
      # Format distribution names
      `Distribution x_1` = dplyr::case_when(
        dist_x1_name == "bernoulli0.5" ~ "Bernoulli(0.5)",
        dist_x1_name == "normal1" ~ "N(0,1)",
        dist_x1_name == "normal5" ~ "N(0,5)",
        TRUE ~ dist_x1_name
      ),
      `Distribution x_2` = dplyr::case_when(
        dist_x2_name == "bernoulli0.5" ~ "Bernoulli(0.5)",
        dist_x2_name == "normal1" ~ "N(0,1)",
        dist_x2_name == "normal5" ~ "N(0,5)",
        TRUE ~ dist_x2_name
      ),

      # Create DGP, Model, and Model Type columns
      DGP = dplyr::case_when(
        problem == "sim_stratified" ~ "STSS DGP",
        problem == "sim_timeScales" ~ "MTS DGP"
      ),
      Model = dplyr::case_when(
        stringr::str_starts(formula_name, "stratified") ~ "STSS PAM",
        stringr::str_starts(formula_name, "timeScales") ~ "MTS PAM"
      ),
      `Model Type` = dplyr::if_else(stringr::str_ends(formula_name, "_ieb"), "$x_1$-only Model", "Full Model"),
      `Model Type` = factor(`Model Type`, levels = c("$x_1$-only Model", "Full Model")),


      # Format other descriptive columns, now including all three transitions
      `Effect Size Scenario` = stringr::str_to_title(beta_scenario),
      Transition = dplyr::case_when(
          transition == "onset"           ~ "Onset",
          transition == "progression_int" ~ "Progression Diff.",
          transition == "progression"     ~ "Progression",
          TRUE                            ~ transition
      )
    ) %>%

    # 3. Rename summary columns
    dplyr::rename(
      `Coefficient (true)` = coef_true,
      `Coefficient (mean)` = coef_mean,
      `Coefficient (lower)` = coef_lower,
      `Coefficient (upper)` = coef_upper,
      `Bias (mean)` = bias_mean,
      `Bias (lower)` = bias_lower,
      `Bias (upper)` = bias_upper,
      `Coverage (mean)` = coverage_mean,
      `Coverage (lower)` = coverage_lower,
      `Coverage (upper)` = coverage_upper
    ) %>%

    # 4. Select and reorder columns for the final CSV file
    dplyr::select(
      `Distribution x_1`,
      `Distribution x_2`,
      `DGP`,
      `Model`,
      `Effect Size Scenario`,
      `Transition`,
      `Model Type`,
      `Coefficient (true)`,
      `Coefficient (mean)`,
      `Coefficient (lower)`,
      `Coefficient (upper)`,
      `Bias (mean)`,
      `Bias (lower)`,
      `Bias (upper)`,
      `Coverage (mean)`,
      `Coverage (lower)`,
      `Coverage (upper)`
    ) %>%

    # 5. Sort the final data frame according to the specified hierarchy
    dplyr::arrange(
      `Distribution x_1`,
      `Distribution x_2`,
      DGP,
      Model,
      `Effect Size Scenario`,
      Transition,
      `Model Type`
    )

  return(csv_output)
}
