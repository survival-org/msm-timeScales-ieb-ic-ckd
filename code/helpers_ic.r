library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(mvtnorm)
library(rlang)
library(survival)
# library(icenReg)
library(mgcv)
library(pammtools)
library(flexsurv)
library(ggplot2)

# data prep ----
wrapper_sim_pexp <- function(
  data,
  job,
  formula = NULL,
  n = 500,
  time_grid = seq(0, 10, by = 0.05),
  ic = TRUE,
  visits_min = 1,
  visits_max = 10,
  ic_mechanism = c("beta", "uniform", "equidistant"),
  round = NULL) {

  df <- tibble::tibble(x1 = rbinom(n, 1, 0.8), x2 = runif(n, 0, 6))

  if(is.null(formula)) {
    formula <- as.formula("~ -3.5 + dgamma(t, 8, 2)*6 - 1.3*x1 + 0.7*x2**2")
  } else {
     formula <- as.formula(formula)
  }

  ndf <- sim_pexp(formula = formula, data = df, cut = time_grid)

  if(ic) {
    ndf <- add_interval_censoring(ndf, max_time = max(ndf$time), visits_min = visits_min, visits_max = visits_max, ic_mechanism = ic_mechanism, round = round) %>%
      filter(time_ic_stop > 0)
  }

  out <- list(ndf, formula)

  return(out)
}

wrapper_sim_weibull <- function(
  data,
  job,
  formula = NULL,
  n = 500,
  time_grid = seq(0, 10, by = 0.05),
  shape = 1.5,
  ic = TRUE,
  visits_min = 1,
  visits_max = 10,
  ic_mechanism = c("beta", "uniform", "equidistant"),
  round = NULL
) {

  # logh(t)=log(α)+(α−1)log(t)−αlog(s) (so base R / English wiki parametrization of Weibull distribution)
  # logh(t)=(time terms) + linear_predictor
  # linear_predictor = -alpha*log(s) --> s = exp(-linear_predictor/alpha)

  # alternatively: use survsim::simple.surv.sim()
  # need to set beta0.cens to a large value to avoid type III right-censoring
  # e.g. formula = "~ -3.5 - 1.3*x1 + sqrt(x2)" and shape = 1.5
  # then anc.ev <- shape; beta0.ev <- -3.5/anc.ev; beta <- c(-1.3/anc.ev, sqrt_x2), where sqrt_x2 has to be created ex ante

  # generate covariates
  df <- tibble::tibble(
    x1 = rbinom(n, 1, 0.8),
    x2 = runif(n, 0, 6)
  )

  # default formula in terms of log hazard
  if (is.null(formula)) {
    formula <- as.formula("~ -3.5 - 1.3*x1 + 0.7*x2**2")
  } else {
    formula <- as.formula(formula)
  }

  # linear predictor (log hazard scale)
  linear_predictor <- with(df, eval(parse(text = as.character(formula)[2])))
  scale <- exp(-linear_predictor / shape)

  # simulate survival times from Weibull distribution
  u <- runif(n)
  survival_times <- scale * (-log(u))^(1 / shape)

  # censoring at maximum follow-up time
  observed_time <- pmin(survival_times, max(time_grid))
  status <- as.integer(survival_times <= max(time_grid))

  # prepare output data frame
  ndf <- tibble::tibble(
    id = 1:n,
    time = observed_time,
    status = status,
    x1 = df$x1,
    x2 = df$x2
  )

  # apply interval censoring if required
  if (ic) {
    ndf <- add_interval_censoring(
      ndf,
      max_time = max(ndf$time),
      visits_min = 1,
      visits_max = 10,
      ic_mechanism = ic_mechanism,
      round = round
    ) %>%
      dplyr::filter(time_ic_stop > 0)
  }

  # generate the symbolic log hazard formula
  log_hazard_formula <- as.formula(
    paste0("~ log(", shape, ") + (", shape, " - 1)*log(t) + ", as.character(formula)[2])
  )

  return(list(ndf = ndf, formula = log_hazard_formula))
}

add_interval_censoring <- function(
  data,
  max_time,
  visits_min = 1,
  visits_max = 10,
  ic_mechanism = c("beta", "uniform", "equidistant"),
  round = NULL) {

  # check necessary columns
  required_cols <- c("id", "time", "status")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: id, time, and status.")
  }
  ic_mechanism <- match.arg(ic_mechanism)

  interval_censor_individual_beta <- function(event_time, event_status, visits_min, visits_max, max_time, round) {
    v <- sample(visits_min:visits_max, 1)
    params <- mvtnorm::rmvnorm(1, mean = c(0, 0), sigma = diag(2))

    shape1 <- abs(params[1]) + 0.5
    shape2 <- abs(params[2]) + 0.5

    quantiles <- seq(0, 1, length.out = v + 2)[-c(1, v + 2)]
    obs_times <- sort(qbeta(quantiles, shape1 = shape1, shape2 = shape2) * max_time)

    if(!is.null(round)) {
      obs_times <- round(obs_times, round)
    }

    obs_times_full <- c(0, obs_times)

    interval_index <- findInterval(event_time, obs_times_full)

    if (interval_index == length(obs_times_full)) {
      # event time exceeds largest observed time: censored after last interval
      time_ic_start <- obs_times_full[interval_index - 1]
      time_ic_stop  <- obs_times_full[interval_index]
      status_ic <- 0
    } else {
      # event time falls between observed interval points
      time_ic_start <- obs_times_full[interval_index]
      time_ic_stop  <- obs_times_full[interval_index + 1]
      status_ic <- ifelse(event_status == 1, 1, 0)
    }

    return(c(time_ic_start = time_ic_start,
            time_ic_stop = time_ic_stop,
            status_ic = status_ic))
  }

  interval_censor_individual_uniform <- function(event_time, event_status, visits_min, visits_max, max_time, round_digits) {
    v <- sample(visits_min:visits_max, 1)

    param <- rnorm(1, mean = 0, sd = 1)
    jitter_sd <- abs(param) + 0.1

    even_times <- seq(0, max_time, length.out = v + 2)[-c(1, v + 2)]
    jittered_times <- even_times + rnorm(v, mean = 0, sd = jitter_sd)

    # Ensure jittered times are strictly within (0, max_time)
    jittered_times <- sort(jittered_times[jittered_times > 0 & jittered_times < max_time])

    # Explicit correction: if no valid jittered times, use max_time as the single interval time
    if (length(jittered_times) == 0) {
      jittered_times <- max_time
    }

    if (!is.null(round_digits)) {
      jittered_times <- round(jittered_times, round_digits)
    }

    obs_times_full <- c(0, jittered_times)

    interval_index <- findInterval(event_time, obs_times_full)

    if (interval_index == length(obs_times_full)) {
      time_ic_start <- obs_times_full[interval_index - 1]
      time_ic_stop  <- obs_times_full[interval_index]
      status_ic <- 0
    } else {
      time_ic_start <- obs_times_full[interval_index]
      time_ic_stop  <- obs_times_full[interval_index + 1]
      status_ic <- ifelse(event_status == 1, 1, 0)
    }

    return(c(time_ic_start = time_ic_start,
            time_ic_stop = time_ic_stop,
            status_ic = status_ic))
  }

  interval_censor_individual_equidistant <- function(event_time, event_status, visits_min, visits_max, max_time, round_digits) {
    breaks <- seq(0, max_time, length.out = visits_max + 1)

    if (!is.null(round_digits)) {
      breaks <- round(breaks, round_digits)
      breaks <- sort(unique(breaks))
    }

    interval_index <- findInterval(event_time, breaks, rightmost.closed = TRUE)

    # Ensure valid interval
    interval_index <- min(max(interval_index, 1), length(breaks) - 1)

    time_ic_start <- breaks[interval_index]
    time_ic_stop  <- breaks[interval_index + 1]
    status_ic <- event_status

    return(c(time_ic_start = time_ic_start,
             time_ic_stop = time_ic_stop,
             status_ic = status_ic))
  }

  # apply to each individual
  interval_fun <- switch(ic_mechanism,
                         beta = interval_censor_individual_beta,
                         uniform = interval_censor_individual_uniform,
                         equidistant = interval_censor_individual_equidistant)

  interval_list <- mapply(interval_fun,
                              event_time = data$time,
                              event_status = data$status,
                              MoreArgs = list(visits_min = visits_min,
                                              visits_max = visits_max,
                                              max_time = max_time,
                                              round = round),
                              SIMPLIFY = FALSE)

  interval_data <- do.call(rbind, interval_list)

  interval_data_df <- as.data.frame(interval_data)

  # return the original data with interval-censored data columns
  final_data <- bind_cols(data, interval_data_df)

  return(final_data)
}

# algorithms ----

## pam wrapper ----
wrapper_pam <- function(
  data,
  job,
  instance,
  bs = "ps",
  k = 10,
  ic_point = c("mid", "end", "mid_end", "exact", "oracle")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  # x_flag <- all(c("x1", "x2") %in% all.vars(formula))
  x_flag <- all(c("x1") %in% all.vars(formula))

  if(ic_point == "mid") {
    df <- df %>%
      mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ .",
    ic_point %in% c("end", "mid_end") ~ "Surv(time_ic_stop, status_ic) ~ .",
    ic_point %in% c("exact", "oracle") ~ "Surv(time, status) ~ .",
    TRUE ~ NA_character_) %>%
    as.formula()

  ped <- as_ped(
    data    = df,
    formula = formula_ped,
    id      = "id")

  if(ic_point == "mid_end") {
    ped <- ped %>%
      mutate(tmid = (tstart + tend) / 2)
  }

  # fit pam
  if (x_flag) {
    if(ic_point == "oracle") {
      formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + x1")
    } else if(ic_point == "mid_end") {
      formula_mod <- paste0("ped_status ~ s(tmid, bs='", bs, "', k=", k, ") + x1")
    } else {
      formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + x1")
    }
  } else {
    if(ic_point == "mid_end") {
      formula_mod <- paste0("ped_status ~ s(tmid, bs='", bs, "', k=", k, ")")
    } else {
      formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ")")
    }
  }

  # if(ic_point != "mid_end") {
  mod <- bam(
    formula = as.formula(formula_mod),
    data    = ped,
    family  = poisson(),
    offset = offset,
    method  = "fREML",
    discrete = TRUE)
  # } else {
  #   mod <- bam(
  #     formula = as.formula(formula_mod),
  #     data    = ped %>% mutate(tend = tmid) %>% select(-tmid),
  #     family  = poisson(),
  #     offset = offset,
  #     method  = "fREML",
  #     discrete = TRUE)
  # }

  # compute coverage for coefficient of x1
  if(x_flag){
    beta_true <- as.numeric(gsub("^.*?([+-]?\\s*\\d+\\.?\\d*)\\s*\\*?\\s*x1.*$", "\\1", gsub("\\s+", "", deparse(formula))))
    summary_mod <- summary(mod)
    beta_est <- summary_mod$p.coeff["x1"]
    beta_se  <- summary_mod$se["x1"]
    ci_lower <- beta_est - 1.96 * beta_se
    ci_upper <- beta_est + 1.96 * beta_se
    beta_cov <- ifelse(beta_true >= ci_lower && beta_true <= ci_upper, 1, 0)

  #   # compute coverage for non-linear spline s(x2)
  #   if (ic_point != "oracle") {

  #     # Define grid of x2-values (from min to max in original data, step size 0.1)
  #     x2_grid <- seq(min(df$x2), max(df$x2), by = 0.1)

  #     # True non-linear effect (sqrt(x2))
  #     sign_coef_func <- strcapture(
  #       "(\\+|\\-)?\\s*(\\d*\\.?\\d*)\\s*\\*?\\s*(\\w+)\\((\\w+)\\)",
  #       paste(deparse(formula), collapse = ""),
  #       proto = list(sign = "+", coef = "1", func = "identity", var = "x2")) |>
  #       transform(
  #         sign = ifelse(sign == "-", -1, 1),
  #         coef = ifelse(coef == "", 1, as.numeric(coef)))

  #     # Define grid of x2-values
  #     x2_grid <- seq(min(df$x2), max(df$x2), by = 0.1)

  #     # Dynamically compute the true nonlinear effect
  #     true_nonlinear_effect <- sign * coef * do.call(func, list(x2_grid))

  #     # Prepare prediction data frame (with median x1 and representative tend)
  #     pred_df <- data.frame(
  #       tend = median(ped$tend),
  #       x1   = median(df$x1),
  #       x2   = x2_grid,
  #       offset = 0
  #     )

  #     # Predict spline term (with standard errors)
  #     spline_pred <- predict(mod, newdata = pred_df, type = "terms", se.fit = TRUE)

  #     # Extract spline estimate and SE
  #     spline_fit <- spline_pred$fit[,"s(x2)"]
  #     spline_se  <- spline_pred$se.fit[,"s(x2)"]

  #     # Compute 95% confidence interval for spline estimate
  #     spline_lower <- spline_fit - 1.96 * spline_se
  #     spline_upper <- spline_fit + 1.96 * spline_se

  #     # Evaluate spline coverage (proportion covered by spline intervals)
  #     spline_cov <-
  #       (true_nonlinear_effect >= spline_lower) &
  #       (true_nonlinear_effect <= spline_upper)

  # ggplot(test_df) +
  #   geom_line(aes(x2_grid, true_nonlinear_effect)) +
  #   geom_line(aes(x2_grid, spline_fit)) +
  #   geom_line(aes(x2_grid, spline_lower), linetype = "dashed") +
  #   geom_line(aes(x2_grid, spline_upper), linetype = "dashed")


  #   } else {

  #     # Oracle scenario: spline coverage is NA
  #     spline_cov <- NA

  #   }
  } else {
    beta_true <- NA
    beta_est <- NA
    beta_se  <- NA
    beta_cov <- NA
  }

  # create new data set
  if(x_flag){
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)),
      x1 = median(df$x1),
      x2 = mean(df$x2)) %>%
      mutate(tmid = (tstart + tend) / 2)
  } else {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend))) %>%
      mutate(tmid = (tstart + tend) / 2)
  }

  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- nd %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)) %>%
    add_hazard(mod, se_mult = qnorm(0.975), ci_type = "default", type = "link") %>%
    rename(loghazard = hazard, loghazard_lower = ci_lower, loghazard_upper = ci_upper) %>%
    add_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    rename(hazard_lower = ci_lower, hazard_upper = ci_upper) %>%
    add_cumu_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    rename(cumu = cumu_hazard) %>%
    add_surv_prob(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    rename(surv = surv_prob) %>%
    mutate(
      loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
      hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
      cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
      surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))

  nd <- nd %>%
    mutate(
      beta_true = beta_true,
      beta_est = beta_est,
      beta_se = beta_se,
      beta_cov = beta_cov) %>%
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov,
      beta_true, beta_est, beta_se, beta_cov)

  return(nd)
}

## cox wrapper ----
wrapper_cox <- function(
  data,
  job,
  instance,
  ic_point = c("mid", "end", "exact", "oracle")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  # x_flag <- all(c("x1", "x2") %in% all.vars(formula))
  x_flag <- all(c("x1") %in% all.vars(formula))

  if(ic_point == "mid") {
    df <- df %>%
      mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ .",
    ic_point == "end" ~ "Surv(time_ic_stop, status_ic) ~ .",
    ic_point %in% c("exact", "oracle") ~ "Surv(time, status) ~ .",
    TRUE ~ NA_character_) %>%
    as.formula()

  ped <- as_ped(
    data    = df,
    formula = formula_ped,
    id      = "id")

  # fit regular coxph model
  if(!x_flag) {
    formula_mod <- update(formula_ped, . ~ 1)
  } else if(ic_point == "oracle") {
    formula_mod <- update(formula_ped, . ~ x1)
  } else {
    formula_mod <- update(formula_ped, . ~ x1)
  }
  mod <- coxph(formula = formula_mod, data = df)

  # compute coverage for coefficient of x1
  if(x_flag) {
    beta_true <- as.numeric(gsub("^.*?([+-]?\\s*\\d+\\.?\\d*)\\s*\\*?\\s*x1.*$", "\\1", gsub("\\s+", "", deparse(formula))))
    summary_mod <- summary(mod)
    beta_est <- summary_mod$coefficients["x1", "coef"]
    beta_se <- coef(summary_mod)["x1", "se(coef)"]
    ci_lower <- beta_est - 1.96 * beta_se
    ci_upper <- beta_est + 1.96 * beta_se
    beta_cov <- ifelse(beta_true >= ci_lower && beta_true <= ci_upper, 1, 0)
  } else {
    beta_true <- NA
    beta_est <- NA
    beta_se  <- NA
    beta_cov <- NA
  }

  # create new dataset exactly as specified
  if(x_flag){
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)),
      x1 = median(df$x1),
      x2 = mean(df$x2)
    )
  } else {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend))
    )
  }

  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- nd %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)
    )

  # get baseline estimates and CI via survfit
  if(x_flag) {
    newdata_pred <- data.frame(x1=median(df$x1), x2=mean(df$x2))
  } else {
    newdata_pred <- data.frame(x1=0, x2=0)
  }
  sf <- survfit(mod, newdata = newdata_pred) # equivalent to survfit(mod)

  # Extract baseline values at matching times
  surv <- approx(sf$time, sf$surv, nd$tend, rule=2)$y
  surv_upper <- approx(sf$time, sf$upper, nd$tend, rule=2)$y
  surv_lower <- approx(sf$time, sf$lower, nd$tend, rule=2)$y

  cumu <- -log(surv)
  cumu_upper <- -log(surv_lower)
  cumu_lower <- -log(surv_upper)

  # baseline_hazard <- c(cumu[1]/nd$tend[1], diff(cumu)/diff(nd$tend))

  # # Linear predictor and its SE
  # if(x_flag) {
  #   lp_pred <- predict(mod, newdata = newdata_pred, type="lp", se.fit=TRUE)
  #   lp <- lp_pred$fit
  #   lp_se <- lp_pred$se.fit
  # }

  # hazard, cumu hazard, and survival adjusted by linear predictor
  nd <- nd %>%
    mutate(
      loghazard = NA_real_,
      hazard = NA_real_, # baseline_hazard * exp(lp): no uncertainty in baseline_hazard!!!
      cumu = cumu,
      surv = surv,
      loghazard_lower = NA_real_,
      loghazard_upper = NA_real_,
      hazard_upper = NA_real_,
      hazard_lower = NA_real_,
      cumu_upper = cumu_upper,
      cumu_lower = cumu_lower,
      surv_upper = surv_upper,
      surv_lower = surv_lower
    ) %>%
    mutate(
      loghazard_cov = NA_real_,
      hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
      cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
      surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))

  nd <- nd %>%
    mutate(
      beta_true = beta_true,
      beta_est = beta_est,
      beta_se = beta_se,
      beta_cov = beta_cov) %>%
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov,
      beta_true, beta_est, beta_se, beta_cov)

  return(nd)
}

## weibull wrapper ----
wrapper_weibull <- function(
  data,
  job,
  instance,
  ic_point = c("mid", "end", "exact", "oracle", "adjustment")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  # x_flag <- all(c("x1", "x2") %in% all.vars(formula))
  x_flag <- all(c("x1") %in% all.vars(formula))

  if(ic_point == "mid") {
    df <- df %>% mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ .",
    ic_point %in% c("end", "adjustment") ~ "Surv(time_ic_stop, status_ic) ~ .",
    ic_point %in% c("exact", "oracle") ~ "Surv(time, status) ~ .",
    TRUE ~ NA_character_
  ) %>% as.formula()

  ped <- as_ped(
    data = df,
    formula = formula_ped,
    id = "id")

  if(ic_point == "adjustment") {
    df_mod <- df %>%
      mutate(
        # Ensure a minimal positive start if zero
        time_ic_start = ifelse(time_ic_start == 0, 1e-6, time_ic_start),
        # For right-censored individuals (status == 0), use the last observed time as the start
        time_ic_start = ifelse(status == 0, time_ic_stop, time_ic_start),
        # And set the upper limit to Inf
        time_ic_stop  = ifelse(status == 0, Inf, time_ic_stop)
      )
  } else {
    df_mod <- df
  }

  if(x_flag) {
    if(ic_point == "adjustment") {
      # Interval-censored formula with covariates
      formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ x1
    } else if(ic_point == "oracle") {
      formula_mod <- update(formula_ped, . ~ x1)
    } else {
      formula_mod <- update(formula_ped, . ~ x1)
    }
  } else {
    if(ic_point == "adjustment") {
      # Interval-censored formula without covariates
      formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ 1
    } else {
      # Regular formula without covariates
      formula_mod <- update(formula_ped, . ~ 1)
    }
  }

  # Model fitting
  mod <- tryCatch(
    withCallingHandlers(
      {
        flexsurvreg(
          formula = formula_mod,
          data = df_mod%>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
          dist = "weibull",
          method = "BFGS",
          inits = c(1, 1, 1))
      },
      warning = function(w) {
        message("Warning detected in BFGS: Falling back to Nelder-Mead")
        invokeRestart("muffleWarning")  # Prevents warning from propagating further
      }
    ),
    error = function(e) {
      message("Error in BFGS: Falling back to Nelder-Mead")
      flexsurvreg(
        formula = formula_mod,
        data = df_mod%>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
        dist = "weibull",
        method = "Nelder-Mead",
        inits = c(1, 1, 1))
    }
  )

  # compute coverage for coefficient of x1
  if(x_flag){
    beta_true <- as.numeric(gsub("^.*?([+-]?\\s*\\d+\\.?\\d*)\\s*\\*?\\s*x1.*$", "\\1", gsub("\\s+", "", deparse(formula))))
    beta_x1_aft       <- mod$res["x1", "est"]
    se_beta_x1_aft    <- mod$res["x1", "se"]
    alpha             <- mod$res["shape", "est"]
    se_log_alpha      <- mod$res.t["shape", "se"]
    beta_est          <- -alpha * beta_x1_aft
    cov_beta_logalpha <- mod$cov["x1", "shape"]
    gradient          <- c(-alpha, -alpha * beta_x1_aft)
    cov_submatrix <- matrix(
      c(se_beta_x1_aft^2, cov_beta_logalpha,
        cov_beta_logalpha, se_log_alpha^2),
      nrow = 2
    )
    beta_var <- t(gradient) %*% cov_submatrix %*% gradient
    beta_se  <- sqrt(beta_var[1, 1])
    ci_lower <- beta_est - 1.96 * beta_se
    ci_upper <- beta_est + 1.96 * beta_se
    beta_cov <- ifelse(beta_true >= ci_lower && beta_true <= ci_upper, 1, 0)
  } else {
    beta_true <- NA
    beta_est <- NA
    beta_se  <- NA
    beta_cov <- NA
  }

  # new data preparation
  if(x_flag) {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)),
      x1 = median(df$x1),
      x2 = mean(df$x2))
  } else {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)))
  }

  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- nd %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)
    )

  # predictions from parametric model
  if(x_flag) {
    newdata_pred <- data.frame(x1=median(df$x1), x2=mean(df$x2))
    pred <- summary(mod, newdata = newdata_pred,
                    type = "hazard", t = nd$tend)
  } else {
    newdata_pred <- data.frame(x1=0, x2=0)
    pred <- summary(mod, newdata = newdata_pred,
                    type = "hazard", t = nd$tend)
  }

  nd <- nd %>%
    mutate(
      loghazard = log(pred[[1]]$est),
      loghazard_lower = log(pred[[1]]$lcl),
      loghazard_upper = log(pred[[1]]$ucl),
      hazard = pred[[1]]$est,
      hazard_lower = pred[[1]]$lcl,
      hazard_upper = pred[[1]]$ucl,
      cumu = cumsum(hazard * intlen),
      cumu_lower = cumsum(hazard_lower * intlen),
      cumu_upper = cumsum(hazard_upper * intlen),
      surv = exp(-cumu),
      surv_lower = exp(-cumu_upper),
      surv_upper = exp(-cumu_lower)
    ) %>%
    mutate(
      loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
      hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
      cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
      surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))

  nd <- nd %>%
    mutate(
      beta_true = beta_true,
      beta_est = beta_est,
      beta_se = beta_se,
      beta_cov = beta_cov) %>%
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov,
      beta_true, beta_est, beta_se, beta_cov)

  return(nd)
}

## generalizedGamma wrapper ----
wrapper_generalizedGamma <- function(
  data,
  job,
  instance,
  ic_point = c("mid", "end", "exact", "oracle", "adjustment")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  # x_flag <- all(c("x1", "x2") %in% all.vars(formula))
  x_flag <- all(c("x1") %in% all.vars(formula))

  if(ic_point == "mid") {
    df <- df %>% mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ .",
    ic_point %in% c("end", "adjustment") ~ "Surv(time_ic_stop, status_ic) ~ .",
    ic_point %in% c("exact", "oracle") ~ "Surv(time, status) ~ .",
    TRUE ~ NA_character_
  ) %>% as.formula()

  ped <- as_ped(
    data = df,
    formula = formula_ped,
    id = "id")

  if(ic_point == "adjustment") {
    df_mod <- df %>%
      mutate(
        # Ensure a minimal positive start if zero
        time_ic_start = ifelse(time_ic_start == 0, 1e-6, time_ic_start),
        # For right-censored individuals (status == 0), use the last observed time as the start
        time_ic_start = ifelse(status == 0, time_ic_stop, time_ic_start),
        # And set the upper limit to Inf
        time_ic_stop  = ifelse(status == 0, Inf, time_ic_stop)
      )
  } else {
    df_mod <- df
  }

  if(x_flag) {
    if(ic_point == "adjustment") {
      # Interval-censored formula with covariates
      formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ x1
    } else if(ic_point == "oracle") {
      formula_mod <- update(formula_ped, . ~ x1)
    } else {
      formula_mod <- update(formula_ped, . ~ x1)
    }
  } else {
    if(ic_point == "adjustment") {
      # Interval-censored formula without covariates
      formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ 1
    } else {
      # Regular formula without covariates
      formula_mod <- update(formula_ped, . ~ 1)
    }
  }

  # Model fitting
  mod <- tryCatch(
    withCallingHandlers(
      {
        flexsurvreg(
          formula = formula_mod,
          data = df_mod%>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
          dist = "gengamma",
          method = "BFGS",
          inits = c(1, 1, 1))
      },
      warning = function(w) {
        message("Warning detected in BFGS: Falling back to Nelder-Mead")
        invokeRestart("muffleWarning")  # Prevents warning from propagating further
      }
    ),
    error = function(e) {
      message("Error in BFGS: Falling back to Nelder-Mead")
      flexsurvreg(
        formula = formula_mod,
        data = df_mod%>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
        dist = "gengamma",
        method = "Nelder-Mead",
        inits = c(1, 1, 1))
    }
  )

  # compute coverage for coefficient of x1
  beta_true <- NA
  beta_est <- NA
  beta_se  <- NA
  beta_cov <- NA

  # new data preparation
  if(x_flag) {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)),
      x1 = median(df$x1),
      x2 = mean(df$x2))
  } else {
    nd <- make_newdata(
      ped,
      tend = sort(unique(ped$tend)))
  }

  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- nd %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)
    )

  # predictions from parametric model
  if(x_flag) {
    newdata_pred <- data.frame(x1=median(df$x1), x2=mean(df$x2))
    pred <- summary(mod, newdata = newdata_pred,
                    type = "hazard", t = nd$tend)
  } else {
    newdata_pred <- data.frame(x1=0, x2=0)
    pred <- summary(mod, type = "hazard", t = nd$tend)
  }

  nd <- nd %>%
    mutate(
      loghazard = log(pred[[1]]$est),
      loghazard_lower = log(pred[[1]]$lcl),
      loghazard_upper = log(pred[[1]]$ucl),
      hazard = pred[[1]]$est,
      hazard_lower = pred[[1]]$lcl,
      hazard_upper = pred[[1]]$ucl,
      cumu = cumsum(hazard * intlen),
      cumu_lower = cumsum(hazard_lower * intlen),
      cumu_upper = cumsum(hazard_upper * intlen),
      surv = exp(-cumu),
      surv_lower = exp(-cumu_upper),
      surv_upper = exp(-cumu_lower)
    ) %>%
    mutate(
      loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
      hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
      cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
      surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))

  nd <- nd %>%
    mutate(
      beta_true = beta_true,
      beta_est = beta_est,
      beta_se = beta_se,
      beta_cov = beta_cov) %>%
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov,
      beta_true, beta_est, beta_se, beta_cov)

  return(nd)
}

# wrapper_generalizedGamma <- function(
#   data,
#   job,
#   instance,
#   ic_point = c("mid", "end", "exact", "oracle", "adjustment")) {

#   ic_point <- match.arg(ic_point)

#   df <- instance[[1]]
#   formula <- instance[[2]]

#   if(ic_point == "mid") {
#     df <- df %>% mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
#   }

#   formula_ped <- case_when(
#     ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
#     ic_point %in% c("end", "adjustment") ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
#     ic_point %in% c("exact", "oracle") ~ "Surv(time, status) ~ x1 + x2",
#     TRUE ~ NA_character_
#   ) %>% as.formula()

#   ped <- as_ped(
#     data = df,
#     formula = formula_ped,
#     id = "id")

#   if(ic_point != "adjustment") {
#     # regular formula
#     formula_mod <- update(formula_ped, . ~ x1 + pspline(x2))
#   }  else {
#     # interval-censored formula
#     formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ x1 + pspline(x2)
#   }

#   # Model fitting
#   mod <- tryCatch(
#     withCallingHandlers(
#       {
#         flexsurvreg(formula_mod, data = df, dist = "gengamma",
#                     inits = c(mu = log(median(df$time)), sigma = 0.5, Q = 1),
#                     method = "BFGS")
#       },
#       warning = function(w) {
#         message("Warning detected in BFGS: Falling back to Nelder-Mead")
#         invokeRestart("muffleWarning")  # Prevents warning from propagating further
#       }
#     ),
#     error = function(e) {
#       message("Error in BFGS: Falling back to Nelder-Mead")
#       flexsurvreg(formula_mod, data = df, dist = "gengamma",
#                   inits = c(mu = log(median(df$time)), sigma = 0.5, Q = 1),
#                   method = "Nelder-Mead")
#     }
#   )

#   # compute coverage for coefficient of x1
#   beta_true <- as.numeric(gsub("^.*?([+-]?\\s*\\d+\\.?\\d*)\\s*\\*?\\s*x1.*$", "\\1", gsub("\\s+", "", deparse(formula))))
#   summary_mod <- summary(mod)
#   beta_est <- 999 # placeholder
#   beta_se <- 999
#   ci_lower <- 999
#   ci_upper <- 999
#   beta_cov <- ifelse(beta_true >= ci_lower && beta_true <= ci_upper, 1, 0)

#   # new data preparation
#   formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

#   nd <- make_newdata(
#     ped,
#     tend = sort(unique(ped$tend)),
#     x1 = median(df$x1),
#     x2 = mean(df$x2)) %>%
#     mutate(
#       loghazard_true = eval(parse(text = formula), envir = pick(everything())),
#       hazard_true = exp(loghazard_true),
#       cumu_true = cumsum(intlen * hazard_true),
#       surv_true = exp(-cumu_true)
#     )

#   # predictions from parametric model
#   pred <- summary(mod, newdata = data.frame(x1 = median(df$x1), x2 = mean(df$x2)),
#                   type = "hazard", t = nd$tend)

#   nd <- nd %>%
#     mutate(
#       loghazard = log(pred[[1]]$est),
#       loghazard_lower = log(pred[[1]]$lcl),
#       loghazard_upper = log(pred[[1]]$ucl),
#       hazard = pred[[1]]$est,
#       hazard_lower = pred[[1]]$lcl,
#       hazard_upper = pred[[1]]$ucl,
#       cumu = cumsum(hazard * intlen),
#       cumu_lower = cumsum(hazard_lower * intlen),
#       cumu_upper = cumsum(hazard_upper * intlen),
#       surv = exp(-cumu),
#       surv_lower = exp(-cumu_upper),
#       surv_upper = exp(-cumu_lower)
#     ) %>%
#     mutate(
#       loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
#       hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
#       cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
#       surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))

#   nd <- nd %>%
#     mutate(
#       beta_true = beta_true,
#       beta_est = beta_est,
#       beta_se = beta_se,
#       beta_cov = beta_cov) %>%
#     select(
#       tstart, tend,
#       loghazard, hazard, cumu, surv,
#       loghazard_true, hazard_true, cumu_true, surv_true,
#       loghazard_cov, hazard_cov, cumu_cov, surv_cov,
#       beta_true, beta_est, beta_se, beta_cov)

#   return(nd)
# }

# postprocessing ----
calc_coverage <- function(data, grouping_vars = NULL, rounding = 3) {

  required_cols <- c(
    grouping_vars,
    "job.id",
    "loghazard_cov", "hazard_cov", "cumu_cov", "surv_cov"
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  jobLevel_results <- data %>%
    {if(!is.null(grouping_vars))
        group_by(., job.id, across(all_of(grouping_vars)))
    else
        group_by(., job.id)} %>%
    summarise(
      loghazard_cov = mean(loghazard_cov, na.rm = TRUE),
      hazard_cov = mean(hazard_cov, na.rm = TRUE),
      cumu_cov = mean(cumu_cov, na.rm = TRUE),
      surv_cov = mean(surv_cov, na.rm = TRUE)
    )

  results <- jobLevel_results %>%
    {if (!is.null(grouping_vars)) group_by(., across(all_of(grouping_vars))) else . } %>%
    summarise(
      `coverage loghazard` = round(mean(loghazard_cov, na.rm = TRUE), rounding),
      `coverage hazard` = round(mean(hazard_cov, na.rm = TRUE), rounding),
      `coverage cumulative hazard` = round(mean(cumu_cov, na.rm = TRUE), rounding),
      `coverage survival probability` = round(mean(surv_cov, na.rm = TRUE), rounding),
      .groups = "drop"
    )

  return(results)
}


calc_rmse <- function(data, grouping_vars = NULL, rounding = 3) {

  required_cols <- c(
    grouping_vars,
    "job.id",
    "loghazard", "hazard", "cumu", "surv",
    "loghazard_true", "hazard_true", "cumu_true", "surv_true"
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  jobLevel_results <- data %>%
    mutate(
      loghazard_sq_error = (loghazard_true - loghazard)**2,
      hazard_sq_error = (hazard_true - hazard)**2,
      cumu_sq_error = (cumu_true - cumu)**2,
      surv_sq_error = (surv_true - surv)**2) %>%
    {if(!is.null(grouping_vars))
        group_by(., job.id, across(all_of(grouping_vars)))
    else
        group_by(., job.id)} %>%
    summarise(
      loghazard_rmse = sqrt(mean(loghazard_sq_error)),
      hazard_rmse = sqrt(mean(hazard_sq_error)),
      cumu_rmse = sqrt(mean(cumu_sq_error)),
      surv_rmse = sqrt(mean(surv_sq_error))
    )

  results <- jobLevel_results %>%
    { if (!is.null(grouping_vars)) group_by(., across(all_of(grouping_vars))) else . } %>%
    summarise(
      `RMSE loghazard` = mean(loghazard_rmse),
      `RMSE hazard` = mean(hazard_rmse),
      `RMSE cumulative hazard` = mean(cumu_rmse),
      `RMSE survival probability` = mean(surv_rmse),
      .groups = "drop"
    )

  return(results)
}


create_linePlot <- function(data, grouping_vars = NULL) {

  required_cols <- c(
    grouping_vars,
    "job.id",
    "tend", "loghazard", "loghazard_true"
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  # Filter out rows with NA loghazard
  plot_data <- data %>%
    filter(!is.na(loghazard))

  # If grouping variables are present, unite them into one 'grouping' col
  if (!is.null(grouping_vars)) {
    plot_data <- plot_data %>%
      unite("grouping", all_of(grouping_vars), remove = FALSE)
  } else {
    plot_data$grouping <- "All"
  }

  unique_groups <- unique(plot_data$grouping)

  plots <- lapply(unique_groups, function(grp) {
    df_grp <- plot_data %>% filter(grouping == grp)

    # Decide line type based on group name
    chosen_line_type <- ifelse(grepl("sim_weibull_weibull|sim_weibull_generalizedGamma|sim_pexp_weibull|sim_pexp_generalizedGamma", grp),
                   "smooth", "step")

    # Base ggplot
    p <- ggplot(df_grp, aes(x = tend, y = loghazard))

    # Draw lines for each job.id according to chosen_line_type
    if (chosen_line_type == "step") {
      p <- p + geom_step(aes(group = job.id), alpha = 0.3)
    } else {
      p <- p + geom_line(aes(group = job.id), alpha = 0.3)
    }

    # Add the true loghazard curve and the GAM-based average estimate
    p +
      geom_line(aes(y = loghazard_true, col = "truth"), lwd = 1.5) +
      geom_smooth(aes(col = "average estimate"), method = "gam",
                  formula = y ~ s(x), se = FALSE) +
      scale_color_brewer("", palette = "Dark2") +
      xlab("time") +
      ylab("loghazard") +
      ggtitle(paste("Group:", grp)) +
      theme_bw()
  })

  names(plots) <- unique_groups
  return(plots)
}


calc_coverage_beta <- function(data, grouping_vars = NULL, rounding = 3) {

  required_cols <- c(
    grouping_vars,
    "job.id",
    "beta_true", "beta_est", "beta_cov"
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  jobLevel_results <- data %>%
    {if(!is.null(grouping_vars))
        group_by(., job.id, across(all_of(grouping_vars)))
    else
        group_by(., job.id)} %>%
    select(all_of(required_cols)) %>%
    slice_head(n = 1) %>%
    mutate(beta_bias = beta_est - beta_true)

  results <- jobLevel_results %>%
    {if (!is.null(grouping_vars)) group_by(., across(all_of(grouping_vars))) else . } %>%
    summarise(
      `coverage beta` = round(mean(beta_cov, na.rm = TRUE), rounding),
      `bias beta` = round(mean(beta_bias, na.rm = TRUE), rounding),
      .groups = "drop"
    )

  return(results)
}

