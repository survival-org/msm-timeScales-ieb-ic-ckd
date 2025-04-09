library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(mgcv)
library(pammtools)
library(flexsurv)
library(ggplot2)

# data prep ----
wrapper_sim <- function(
  data,
  job,
  formulas_list = list(
    list(from = 0, to = 1, formula = ~ -3.5 + dgamma(t, 8, 2)*6 + 1.3*x1 + 0.8*x2),
    list(from = 0, to = "death", formula = ~ -3.0 + dgamma(t, 8, 2)*6 + 1.3*x1 + 0.8*x2),
    list(1, "death", ~ -2.1 + 1.3*x1 + 0.8*x2)
  ),
  terminal_states = c("death"),
  cut = seq(0, 5, by = 0.01),
  n = 2500,
  round = 2,
  cens_type = c("none", "right", "interval", "left"),
  cens_dist = c("weibull", "exponential", "lognormal", "uniform"),
  cens_params = c(1, 1)){

  data <- cbind.data.frame(
    id = 1:n,
    x1 = rbinom(n, 1, 0.5),
    x2 = rbinom(n, 1, 0.5),
    x3 = runif(n, 0, 3),
    from = 0,
    t = 0)

  df <- sim_pexp_msm(formulas_list, data, cut, terminal_states, round = round, add_counterfactuals = FALSE) %>%
    mutate(transition = as.factor(transition))

  if(cens_type != "none") {
    df <- df %>%
      add_cens(type = cens_type, distribution = cens_dist, parameters = cens_params, round = round)
  }

  ped <- as_ped_multistate(
    data       = df %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

  out <- list(
    df = df,
    ped = ped
  )

  return(out)
}

wrapper_msm <- function(
  data,
  job,
  instance,
  formula = "ped_status ~ s(tend, by=transition) + transition + x1*transition") {

  ped <- instance$ped

  mod <- bam(as.formula(formula)
                      , data = ped
                      , family=poisson()
                      , offset=offset
                      , discrete = T
                      , method = "fREML"
  )
  summary_mod <- summary(mod)

  # Create a dataframe with coefficient information
  coef_df <- data.frame(
    transition = c("0->1 (onset)", "1->death (progression)"),
    coefficient = c(summary_mod$p.coeff["x1"],
                   summary_mod$p.coeff["transition1->death:x1"]),
    std_error = c(summary_mod$se["x1"],
                 summary_mod$se["transition1->death:x1"]),
    p_value = c(summary_mod$p.pv["x1"],
               summary_mod$p.pv["transition1->death:x1"]),
    stringsAsFactors = FALSE
  )

  # Round numeric values for better readability
  coef_df$coefficient <- round(coef_df$coefficient, 4)
  coef_df$std_error <- round(coef_df$std_error, 4)
  coef_df$p_value <- round(coef_df$p_value, 4)

  return(coef_df)
}

wrapper_cor <- function(
  data,
  job,
  instance,
  vars = c("x1", "x2", "x3")) {

  df <- instance$df

  # Check if all variables are in the dataframe
  if (!all(vars %in% colnames(df))) {
    missing_vars <- vars[!vars %in% colnames(df)]
    stop("The following variables are not found in the data: ",
         paste(missing_vars, collapse = ", "))
  }

  cor_df <- calc_statewise_correlation(
    data = df,
    from_states = unique(df$from),
    vars = c("x1", "x2", "x3")
  )

  return(cor_df)
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
