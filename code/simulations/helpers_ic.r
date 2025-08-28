library(data.table)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(stringr)
library(mvtnorm)
library(rlang)
library(survival)
library(mgcv)
library(pammtools)
library(flexsurv)
library(ggplot2)
theme_set(theme_bw())
library(patchwork)

# data prep ----
wrapper_sim_pexp <- function(
  data,
  job,
  formula = NULL,
  n = 500,
  cut = seq(0, 10, by = 0.1),
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

  ndf <- sim_pexp(formula = formula, data = df, cut = cut)

  if(ic) {
    ndf <- add_interval_censoring(ndf, max_time = max(ndf$time), visits_min = visits_min, visits_max = visits_max, ic_mechanism = ic_mechanism, round = round) %>%
      filter(time_ic_stop > 0)
  } else {
    ndf <- ndf %>%
      mutate(time = round(time, round))
  }

  out <- list(ndf, formula)

  return(out)
}


wrapper_sim_weibull <- function(
  data,
  job,
  formula = NULL,
  n = 500,
  cut = seq(0, 10, by = 0.1),
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
  observed_time <- pmin(survival_times, max(cut))
  status <- as.integer(survival_times <= max(cut))

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
      visits_min = visits_min,
      visits_max = visits_max,
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


wrapper_sim_icenReg <- function(
  data,
  job,
  n = 500,
  cut = seq(0, 10, by = 0.1),
  scale = 2,
  shape = 1.5,
  inspections = 20,
  inspectLength = 4,
  round = NULL) {

  max_time <- max(cut)

  ndf <- simIC_weib(
    n = n,
    b1 = 0,
    b2 = 0,
    model = "ph",
    shape = shape,
    scale = scale,
    inspections = inspections,
    inspectLength = inspectLength,
    rndDigits = round,
    prob_cen = 1
  ) %>%
    filter(l <= max_time) %>%
    mutate(
      u = ifelse(u == Inf, l, u),
      u = pmin(u, max_time),
      status = ifelse(time > max_time | l == u, 0, 1),
      status_ic = ifelse(time > max_time | l == u, 0, 1),
      time = ifelse(l == u, u, time),
      time = pmin(time, max_time),
    ) %>%
    rename(
      time_ic_start = l,
      time_ic_stop = u
    ) %>%
    filter(time_ic_stop > time_ic_start) # otherwise incompatible with model algorithms if length(interval)=0

  formula_call <- substitute(
    ~ log(alpha) - log(sigma) + (alpha - 1) * (log(t) - log(sigma)),
    list(alpha = shape, sigma = scale)
  )
  formula <- as.formula(formula_call)

  out <- list(ndf, formula)

  return(out)
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

simIC_weib <- function (n = 100, b1 = 0.5, b2 = -0.5, model = "ph", shape = 2,
    scale = 2, inspections = 2, inspectLength = 2.5, rndDigits = NULL,
    prob_cen = 1){
    rawQ <- runif(n)
    x1 <- runif(n, -1, 1)
    x2 <- 1 - 2 * rbinom(n, 1, 0.5)
    nu <- exp(x1 * b1 + x2 * b2)
    if (model == "ph")
        adjFun <- function(x, nu) {
            1 - x^(1/nu)
        }
    else if (model == "po")
        adjFun <- function(x, nu) {
            1 - x * (1/nu)/(x * 1/nu - x + 1)
        }
    adjQ <- adjFun(rawQ, nu)
    trueTimes <- qweibull(adjQ, shape = shape, scale = scale)
    obsTimes <- runif(n = n, max = inspectLength)
    if (!is.null(rndDigits))
        obsTimes <- round(obsTimes, rndDigits)
    l <- rep(0, n)
    u <- rep(0, n)
    caught <- trueTimes < obsTimes
    u[caught] <- obsTimes[caught]
    l[!caught] <- obsTimes[!caught]
    if (inspections > 1) {
        for (i in 2:inspections) {
            oldObsTimes <- obsTimes
            obsTimes <- oldObsTimes + runif(n, max = inspectLength)
            if (!is.null(rndDigits))
                obsTimes <- round(obsTimes, rndDigits)
            caught <- trueTimes >= oldObsTimes & trueTimes <
                obsTimes
            needsCatch <- trueTimes > obsTimes
            u[caught] <- obsTimes[caught]
            l[needsCatch] <- obsTimes[needsCatch]
        }
    }
    else {
        needsCatch <- !caught
    }
    u[needsCatch] <- Inf
    if (sum(l > u) > 0)
        stop("warning: l > u! Bug in code")
    # isCensored <- rbinom(n = n, size = 1, prob = prob_cen) ==
    #     1
    # l[!isCensored] <- trueTimes[!isCensored]
    # u[!isCensored] <- trueTimes[!isCensored]
    if (sum(l == Inf) > 0) {
        allTimes <- c(l, u)
        allFiniteTimes <- allTimes[allTimes < Inf]
        maxFiniteTime <- max(allFiniteTimes)
        l[l == Inf] <- maxFiniteTime
    }
    return(data.frame(l = l, u = u, time = trueTimes, x1 = x1, x2 = x2))
}

# algorithms ----

## pam wrapper ----
wrapper_pam <- function(
  data,
  job,
  instance,
  bs = "ps",
  k = 20,
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
    cut     = data$cut,
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
    cut     = data$cut,
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
  ic_point = c("mid", "end", "exact", "oracle", "adjustment"),
  fct = c("flexsurvreg", "survreg")
  ) {

  ic_point <- match.arg(ic_point)
  fct <- match.arg(fct)

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
    cut     = data$cut,
    id = "id")

  if(ic_point == "adjustment") {
    df_mod <- df %>%
      mutate(
        # Ensure a minimal positive start if zero
        time_ic_start = ifelse(time_ic_start == 0, 1e-6, time_ic_start),
        # For right-censored individuals according to IC scheme, use the last observed time as the start
        time_ic_start = ifelse(status_ic == 0, time_ic_stop, time_ic_start),
        # And set the upper limit to Inf
        time_ic_stop  = ifelse(status_ic == 0, Inf, time_ic_stop)
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
  if(fct == "flexsurvreg") {
    mod <- tryCatch(
      withCallingHandlers(
        {
          flexsurvreg(
            formula = formula_mod,
            data = df_mod %>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
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
          data = df_mod %>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
          dist = "weibull",
          method = "Nelder-Mead",
          inits = c(1, 1, 1))
      }
    )
  } else {
    mod <- survreg(
      formula = formula_mod,
      data = df_mod %>% mutate(time_ic_start = ifelse(time_ic_start==0, 1e-6, time_ic_start)),
      dist = "weibull")
  }

  # compute coverage for coefficient of x1
  if(x_flag){

    beta_true <- as.numeric(gsub("^.*?([+-]?\\s*\\d+\\.?\\d*)\\s*\\*?\\s*x1.*$", "\\1", gsub("\\s+", "", deparse(formula))))

    if(fct == "flexsurvreg") {
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
      stop("survreg not implemented yet for coverage of covariates")
    }

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
    newdata_pred <- data.frame(x1 = median(df$x1), x2 = mean(df$x2))
    if(fct == "flexsurvreg") {
      pred <- summary(mod, newdata = newdata_pred,
                      type = "hazard", t = nd$tend)
    } else {
      # For survreg: compute hazard manually
      # Get linear predictor and sigma
      lp <- predict(mod, newdata = newdata_pred, type = "lp")
      sigma <- mod$scale
      a <- 1/sigma        # shape parameter
      lambda <- exp(lp)   # scale parameter
      # Compute hazard at each time point in nd$tend
      hazard <- (a / lambda) * (nd$tend / lambda)^(a - 1)
      # Note: Without additional delta-method calculations, we do not get SE/CI.
      # For consistency downstream, we create a list with est, lcl, and ucl all equal to hazard.
      pred <- list(list(est = hazard,
                        lcl = hazard,  # placeholder: no CI available for survreg predictions here
                        ucl = hazard))
    }
  } else {
    newdata_pred <- data.frame(x1 = 0, x2 = 0)
    if(fct == "flexsurvreg") {
      pred <- summary(mod, newdata = newdata_pred,
                      type = "hazard", t = nd$tend)
    } else {
      lp <- predict(mod, newdata = newdata_pred, type = "lp")
      sigma <- mod$scale
      a <- 1/sigma
      lambda <- exp(lp)
      hazard <- (a / lambda) * (nd$tend / lambda)^(a - 1)
      pred <- list(list(est = hazard, lcl = hazard, ucl = hazard))
    }
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
    cut     = data$cut,
    id = "id")

  if(ic_point == "adjustment") {
    df_mod <- df %>%
      mutate(
        # Ensure a minimal positive start if zero
        time_ic_start = ifelse(time_ic_start == 0, 1e-6, time_ic_start),
        # For right-censored individuals according to IC scheme, use the last observed time as the start
        time_ic_start = ifelse(status_ic == 0, time_ic_stop, time_ic_start),
        # And set the upper limit to Inf
        time_ic_stop  = ifelse(status_ic == 0, Inf, time_ic_stop)
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
calc_pointwise_coverage_bh <- function(data,
                                    grouping_vars = NULL,
                                    cut_var = "tend",
                                    rounding = 3) {
  # 1) check
  req <- c(grouping_vars, "job.id", cut_var,
           "loghazard_cov", "hazard_cov",
           "cumu_cov", "surv_cov")
  if (!all(req %in% names(data))) {
    stop("Data is missing required columns.")
  }

  # 2) group by group + cut
  grp_cols <- c(grouping_vars, cut_var)
  results <- data %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(
      coverage_loghazard = mean(loghazard_cov, na.rm = TRUE),
      ci_loghazard = list(
        binom.test(
          x = sum(loghazard_cov == 1, na.rm = TRUE),
          n = n(),
          p = 0.5
        )$conf.int
      ),
      coverage_hazard    = mean(hazard_cov, na.rm = TRUE),
      ci_hazard = list(
        binom.test(
          x = sum(hazard_cov == 1, na.rm = TRUE),
          n = n(),
          p = 0.5
        )$conf.int
      ),
      coverage_cumu      = mean(cumu_cov, na.rm = TRUE),
      ci_cumu = list(
        binom.test(
          x = sum(cumu_cov == 1, na.rm = TRUE),
          n = n(),
          p = 0.5
        )$conf.int
      ),
      coverage_surv      = mean(surv_cov, na.rm = TRUE),
      ci_surv = list(
        binom.test(
          x = sum(surv_cov == 1, na.rm = TRUE),
          n = n(),
          p = 0.5
        )$conf.int
      ),
      .groups = "drop"
    ) %>%
    # 3) pull out lower/upper and round
    mutate(
      across(starts_with("coverage"), ~ round(.x, rounding)),
      coverage_loghazard_lower = round(map_dbl(ci_loghazard, 1), rounding),
      coverage_loghazard_upper = round(map_dbl(ci_loghazard, 2), rounding),
      coverage_hazard_lower    = round(map_dbl(ci_hazard, 1), rounding),
      coverage_hazard_upper    = round(map_dbl(ci_hazard, 2), rounding),
      coverage_cumu_lower      = round(map_dbl(ci_cumu, 1), rounding),
      coverage_cumu_upper      = round(map_dbl(ci_cumu, 2), rounding),
      coverage_surv_lower      = round(map_dbl(ci_surv, 1), rounding),
      coverage_surv_upper      = round(map_dbl(ci_surv, 2), rounding)
    ) %>%
    select(-ci_loghazard, -ci_hazard, -ci_cumu, -ci_surv)

  return(results)
}


calc_pointwise_rmse_bh <- function(data,
                                    grouping_vars = NULL,
                                    cut_var = "tend",
                                    rounding = 3) {
  # 1) check required columns
  req_cols <- c(grouping_vars,
                "job.id", cut_var,
                "loghazard", "hazard", "cumu", "surv",
                "loghazard_true", "hazard_true", "cumu_true", "surv_true")
  if (!all(req_cols %in% names(data))) {
    stop("Data is missing required columns.")
  }

  # 2) drop infinite predictions, then group + compute RMSE
  grp_cols <- c(grouping_vars, cut_var)
  results <- data %>%
    # drop rows where any of the four predictions is ±Inf
    filter(
      is.finite(loghazard),
      is.finite(hazard),
      is.finite(cumu),
      is.finite(surv)
    ) %>%
    # compute squared errors
    mutate(
      loghazard_se = (loghazard - loghazard_true)^2,
      hazard_se    = (hazard    - hazard_true   )^2,
      cumu_se      = (cumu      - cumu_true     )^2,
      surv_se      = (surv      - surv_true     )^2
    ) %>%
    group_by(across(all_of(grp_cols))) %>%
    summarise(
      rmse_loghazard = sqrt(mean(loghazard_se, na.rm = TRUE)),
      rmse_hazard    = sqrt(mean(hazard_se,    na.rm = TRUE)),
      rmse_cumu      = sqrt(mean(cumu_se,      na.rm = TRUE)),
      rmse_surv      = sqrt(mean(surv_se,      na.rm = TRUE)),
      .groups        = "drop"
    ) %>%
    # 3) round
    mutate(across(starts_with("rmse_"), ~ round(.x, rounding)))

  return(results)
}


create_linePlot <- function(
  data,
  grouping_vars = NULL,
  scale = c("loghazard", "hazard", "cumulativehazard", "survivalfunction"),
  font_size = 14,
  alpha = 1,
  show_title = TRUE  # NEW
) {
  match.arg(scale)

  required_cols <- c(grouping_vars, "job.id", "tend", "loghazard", "loghazard_true")
  if (!all(required_cols %in% colnames(data))) stop("Data is missing required columns.")

  plot_data <- data %>% dplyr::filter(!is.na(loghazard))

  if (!is.null(grouping_vars)) {
    plot_data <- plot_data %>% tidyr::unite("grouping", dplyr::all_of(grouping_vars), remove = FALSE)
  } else {
    plot_data$grouping <- "All"
  }

  unique_groups <- unique(plot_data$grouping)

  plots <- lapply(unique_groups, function(grp) {
    df_grp <- plot_data %>% dplyr::filter(grouping == grp)

    chosen_line_type <- ifelse(grepl("sim_weibull_weibull|sim_weibull_generalizedGamma|sim_pexp_weibull|sim_pexp_generalizedGamma|sim_icenReg_weibull|sim_icenReg_generalizedGamma", grp),
                               "smooth", "step")

    mapping <- list(
      loghazard        = list(var = "loghazard", label = "Log-hazard"),        # CHANGED label
      hazard           = list(var = "hazard",    label = "Hazard"),            # CHANGED label
      cumulativehazard = list(var = "cumu",      label = "Cumulative hazard"), # CHANGED label
      survivalfunction = list(var = "surv",      label = "Survival function")  # CHANGED label
    )
    y_info     <- mapping[[scale]]
    y_var      <- y_info$var
    y_var_true <- paste0(y_var, "_true")
    ylab       <- y_info$label

    rng     <- range(df_grp[[y_var]], df_grp[[y_var_true]], na.rm = TRUE)
    rng_min <- rng[1]
    rng_span<- diff(rng)
    ref_y   <- rng_min + 0.95 * rng_span

    cov_line <- df_grp %>%
      dplyr::group_by(tend) %>%
      dplyr::summarise(
        cov_mean = mean(.data[[paste0(y_var, "_cov")]], na.rm = TRUE),
        .groups  = "drop"
      ) %>%
      dplyr::mutate(cov_y = rng_min + cov_mean * rng_span)

    p <- ggplot2::ggplot(df_grp, ggplot2::aes(x = tend, y = .data[[y_var]]))
    truth <- ggplot2::geom_line(ggplot2::aes(y = .data[[y_var_true]], col = "truth"), linewidth = 1.5)

    if (chosen_line_type == "step") {
      p <- p + ggplot2::geom_step(ggplot2::aes(group = job.id), alpha = alpha)
    } else {
      p <- p + ggplot2::geom_line(ggplot2::aes(group = job.id), alpha = alpha)
    }

    p <- p +
      truth +
      ggplot2::geom_smooth(ggplot2::aes(col = "average estimate"),
                           method = "gam", formula = y ~ s(x), se = FALSE) +
      ggplot2::geom_line(data = cov_line, ggplot2::aes(x = tend, y = cov_y),
                         color = "blue", linetype = "dotted") +
      ggplot2::geom_hline(yintercept = ref_y, linetype = "dashed", color = "red") +
      ggplot2::scale_color_brewer("", palette = "Dark2") +
      ggplot2::scale_y_continuous(
        name = ylab,
        sec.axis = ggplot2::sec_axis(
          trans  = ~ (. - rng_min) / rng_span,
          name   = "Coverage",                                         # CHANGED name
          breaks = seq(0, 0.8, 0.2),                                   # CHANGED ticks
          labels = scales::percent_format(accuracy = 1)                # CHANGED labels
        )
      ) +
      ggplot2::xlab("Time") +                                         # CHANGED x-axis title
      ggplot2::theme_bw() +
      ggplot2::theme(
        legend.position = c(0.05, 0.95),
        legend.justification = c("left", "top"),
        legend.text  = ggplot2::element_text(size = font_size),
        legend.title = ggplot2::element_text(size = font_size),
        axis.text    = ggplot2::element_text(size = font_size),
        axis.title   = ggplot2::element_text(size = font_size)
      )

    if (show_title) {
      p <- p + ggplot2::ggtitle(paste("Group:", grp))                 # keep old default behavior
    }
    p
  })

  names(plots) <- unique_groups
  return(plots)
}



plot_coverage_bh <- function(data,
                             grouping_vars = NULL,
                             scale = c("loghazard", "hazard",
                                       "cumulativehazard", "survivalfunction"),
                             font_size = 14,
                             dot_size = 2.5) {
  scale <- match.arg(scale)

  var_main  <- switch(scale,
                      loghazard          = "coverage_loghazard",
                      hazard             = "coverage_hazard",
                      cumulativehazard   = "coverage_cumu",
                      survivalfunction   = "coverage_surv")
  var_lower <- paste0(var_main, "_lower")
  var_upper <- paste0(var_main, "_upper")

  data <- data %>%
    dplyr::filter(!is.na(.data[[var_main]]))

  groups <- if (!is.null(grouping_vars)) {
    split(data, data[grouping_vars], drop = TRUE)
  } else {
    list(overall = data)
  }

  plots <- lapply(names(groups), function(name) {
    df <- groups[[name]]
    # df$algorithm <- factor(df$algorithm, levels = unique(df$algorithm))
    df$algorithm <- droplevels(df$algorithm)
    pos <- position_dodge(width = 0.6)

    ggplot(df, aes(x = algorithm,
                   y = .data[[var_main]],
                   colour = ic_point)) +
      geom_errorbar(aes(ymin = .data[[var_lower]],
                        ymax = .data[[var_upper]]),
                    position = pos,
                    width = 0.15,
                    linewidth = 1) +     # CHANGED: thickness
      geom_point(position = pos, size = dot_size) +  # CHANGED: dot size
      geom_hline(yintercept = 0.95, colour = "red", linetype = "dashed") +
      scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        labels = scales::percent_format(accuracy = 1)
      ) +
      labs(title = NULL, x = NULL, y = "Coverage", colour = "Estimation Point") +
      theme_bw() +
      theme(
        # axis.text.x  = element_text(size = font_size),
        axis.text.x = element_text(size = font_size, angle = 30, hjust = 1),
        axis.text.y  = element_text(size = font_size),
        axis.title.x = element_text(size = font_size),
        axis.title.y = element_text(size = font_size),
        legend.text  = element_text(size = font_size),
        legend.title = element_text(size = font_size)
      )
  })

  names(plots) <- names(groups)
  plots
}



summarize_fe <- function(data, grouping_vars = NULL, rounding = 3) {
  # check required columns
  req <- c(grouping_vars, "job.id", "beta_true", "beta_est", "beta_se", "beta_cov")
  if (!all(req %in% names(data))) {
    stop("Data is missing required columns.")
  }

  # collapse to one row per job
  jobLevel <- data %>%
    { if (!is.null(grouping_vars))
        group_by(., job.id, across(all_of(grouping_vars)))
      else
        group_by(., job.id)
    } %>%
    slice_head(n = 1) %>%
    mutate(beta_bias = beta_est - beta_true)

  # summarise per group
  res <- jobLevel %>%
    { if (!is.null(grouping_vars))
        group_by(., across(all_of(grouping_vars)))
      else
        .
    } %>%
    summarise(
      coverage   = mean(beta_cov, na.rm = TRUE),
      bias       = mean(beta_bias, na.rm = TRUE),
      beta_est   = mean(beta_est,  na.rm = TRUE),
      beta_se    = mean(beta_se,   na.rm = TRUE),
      ci         = list(
                     binom.test(
                       x = sum(beta_cov, na.rm = TRUE),
                       n = length(beta_cov),
                       p = 0.5
                     )$conf.int
                   ),
      .groups = "drop"
    ) %>%
    mutate(
      coverage       = round(coverage, rounding),
      bias           = round(bias,     rounding),
      beta_est       = round(beta_est, rounding),
      beta_se        = round(beta_se,  rounding),
      coverage_lower = round(map_dbl(ci, 1), rounding),
      coverage_upper = round(map_dbl(ci, 2), rounding)
    ) %>%
    select(-ci)

  return(res)
}

plot_coef_fe <- function(data, grouping_vars = NULL, font_size = 14) {  # new font_size arg
  # 1) split into groups (or overall)
  groups <- if (!is.null(grouping_vars)) {
    split(data, data[grouping_vars], drop = TRUE)
  } else {
    list(overall = data)
  }

  # 2) one grouped boxplot per group
  plots <- lapply(names(groups), function(name) {
    df <- groups[[name]]
    true_val <- unique(df$beta_true)
    if (length(true_val) != 1)
      stop("Group ", name, " has multiple beta_true values.")

    ggplot(df, aes(x = factor(algorithm),
                   y = beta_est,
                   fill = ic_point)) +
      geom_boxplot(position = position_dodge(width = 0.8)) +
      geom_hline(yintercept = true_val,
                 color = "red",
                 linetype = "dashed") +
      labs(
        title = NULL,                  # removed title
        x     = NULL,                  # removed x-axis title
        y     = "Estimated Coefficient",
        fill  = "Estimation Time Point"  # updated legend title
      ) +
      scale_y_continuous(
        breaks = seq(-2.0, -0.5, 0.5),    # fixed y-axis labels
        limits = c(-2.1, -0.45)           # fixed y-axis limits
      ) +
      theme_bw() +
      theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = font_size),
        axis.text.y  = element_text(size = font_size),
        axis.title.y = element_text(size = font_size),
        legend.text  = element_text(size = font_size),
        legend.title = element_text(size = font_size)
      )
  })

  names(plots) <- names(groups)
  plots
}



plot_coverage_fe <- function(data, grouping_vars = NULL) {
  # 1) split into groups (or overall)
  groups <- if (!is.null(grouping_vars)) {
    split(data, data[grouping_vars], drop = TRUE)
  } else {
    list(overall = data)
  }

  # 2) one barplot per group
  plots <- lapply(names(groups), function(name) {
    df <- groups[[name]]
    # ensure algorithm is a factor
    df$algorithm <- factor(df$algorithm, levels = unique(df$algorithm))

    ggplot(df, aes(x = algorithm, y = coverage, fill = ic_point)) +
      geom_col(position = position_dodge(width = 0.8)) +
      geom_errorbar(aes(ymin = coverage_lower, ymax = coverage_upper),
                    position = position_dodge(width = 0.8),
                    width = 0.2) +
      geom_hline(yintercept = 0.95, color = "orange", linetype = "dashed") +
      labs(
        title = name,
        x     = "Algorithm",
        y     = "Coverage",
        fill  = "IC Point"
      ) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  names(plots) <- names(groups)
  plots
}


# tables ----

## bh ----
create_coverage_summary_bh <- function(data, coverage_type, rounding_digits = 2) {

  # --- Step 1: Dynamic column names and format string ---
  cov_col <- paste0("coverage_", coverage_type)
  lower_col <- paste0("coverage_", coverage_type, "_lower")
  upper_col <- paste0("coverage_", coverage_type, "_upper")

  if (!all(c(cov_col, lower_col, upper_col) %in% names(data))) {
    stop(paste("Columns for coverage_type '", coverage_type, "' not found.", sep=""))
  }

  # Create dynamic format string for sprintf based on rounding_digits
  sprintf_format <- paste0("%.", rounding_digits, "f\\n(%.", rounding_digits, "f; %.", rounding_digits, "f)")

  # --- Step 2: Capture factor levels and process data ---
  ic_point_levels <- levels(data$ic_point)
  ic_mechanism_levels <- levels(data$ic_mechanism)

  processed_data <- data %>%
    mutate(
      # Create full DGP names for column headers
      problem_label = case_when(
        problem == "sim_pexp" ~ "Piecewise Exponential DGP",
        problem == "sim_weibull" ~ "Weibull DGP",
        problem == "sim_icenReg" ~ "icenReg DGP",
        TRUE ~ as.character(problem)
      ),
      algorithm = case_when(
        algorithm == "pam" ~ "PAM",
        algorithm == "generalizedGamma" ~ "GenGamma",
        TRUE ~ str_to_title(as.character(algorithm))
      ),
      ic_point = str_to_title(as.character(ic_point)),
      ic_mechanism = str_to_title(as.character(ic_mechanism)),
      # Create the formatted string, robustly handling NA/NaN
      coverage_formatted = if_else(
        is.na(.data[[cov_col]]) | is.na(.data[[lower_col]]) | is.na(.data[[upper_col]]),
        NA_character_,
        sprintf(sprintf_format, .data[[cov_col]], .data[[lower_col]], .data[[upper_col]])
      )
    )

  # --- Step 3: Pivot to create the wide, multi-DGP structure ---
  wide_table <- processed_data %>%
    pivot_wider(
      id_cols = c(algorithm, ic_point),
      names_from = c(problem_label, ic_mechanism),
      values_from = coverage_formatted
    )

  # --- Step 4: Dynamically order columns based on DGP and mechanism factors ---
  dgp_order <- c("Piecewise Exponential DGP", "Weibull DGP", "icenReg DGP")
  mechanism_order <- str_to_title(ic_mechanism_levels)

  # Create all possible column names in the correct final order
  ordered_data_cols <- as.vector(t(outer(dgp_order, mechanism_order, paste, sep = "_")))

  # Filter for only those columns that actually exist in the pivoted data
  final_col_order <- c("algorithm", "ic_point", intersect(ordered_data_cols, names(wide_table)))

  # --- Step 5: Arrange, filter NA rows, and return ---
  final_table <- wide_table %>%
    select(all_of(final_col_order)) %>%
    # New step: Drop rows where all data columns are NA
    filter(if_any(-c(algorithm, ic_point), ~ !is.na(.))) %>%
    arrange(
      factor(algorithm, levels = c("PAM", "Cox", "Weibull", "GenGamma")),
      factor(ic_point, levels = str_to_title(ic_point_levels))
    )

  return(final_table)
}


generate_multi_dgp_latex <- function(summary_df, coverage_type) {

  # (The escape_latex helper function remains the same and is correct)
  escape_latex <- function(text) {
    if (is.na(text)) return("")
    text %>%
      str_replace_all("&", "\\\\&") %>% str_replace_all("%", "\\\\%") %>%
      str_replace_all("_", "\\\\_")
  }

  # --- (Step 1: Dynamic header generation logic remains identical) ---
  caption_line <- paste0("\\caption{\\captionicbhcoverage", coverage_type, "}")
  label_line <- paste0("\\label{tab:sim-ic-bh-coverage-", coverage_type, "}")
  data_cols <- names(summary_df)[-c(1, 2)]
  col_info <- tibble(
    col_name = data_cols,
    dgp = sapply(strsplit(col_name, "_"), `[`, 1),
    mechanism = sapply(strsplit(col_name, "_"), `[`, 2)
  )
  dgp_structure <- col_info %>%
    group_by(dgp) %>%
    summarise(n = n(), .groups = 'drop') %>%
    mutate(dgp = factor(dgp, levels = c("Piecewise Exponential DGP", "Weibull DGP", "Interval Censoring DGP"))) %>%
    arrange(dgp)
  multicolumn_items <- mapply(function(dgp, n) {
    paste0("\\multicolumn{", n, "}{c}{", dgp, "}")
  }, dgp_structure$dgp, dgp_structure$n)
  multicolumn_line <- paste("& &", paste(multicolumn_items, collapse = " & "), "\\\\")
  cmid_start_col <- 3
  cmidrule_items <- c()
  for (n_cols in dgp_structure$n) {
    cmid_end_col <- cmid_start_col + n_cols - 1
    cmidrule_items <- c(cmidrule_items, paste0("\\cmidrule(lr){", cmid_start_col, "-", cmid_end_col, "}"))
    cmid_start_col <- cmid_end_col + 1
  }
  cmidrule_line <- paste(cmidrule_items, collapse = " ")
  mechanism_header_line <- paste(
    "Algorithm & \\shortstack[l]{Estimation \\\\ Time Point}",
    paste(col_info$mechanism, collapse = " & "),
    sep = " & "
  )

  # --- Step 2: Assemble the full table with the key change ---
  num_total_data_cols <- length(data_cols)
  tabular_def <- paste0("\\begin{tabular}{ll", paste(rep("c", num_total_data_cols), collapse = ""), "}")

  header <- c(
    "\\begin{table*}[htbp]", "\\centering",
    caption_line,
    label_line,
    "\\begin{sideways}",
    # *** THIS IS THE KEY CHANGE TO MAKE THE TABLE NARROWER ***
    # We set a smaller column separation. Default is 6pt.
    "\\setlength{\\tabcolsep}{4pt}",
    tabular_def, "\\toprule",
    multicolumn_line, cmidrule_line, paste0(mechanism_header_line, " \\\\"), "\\midrule"
  )

  # --- (Body generation and footer remain identical) ---
  body_rows <- c()
  for (i in 1:nrow(summary_df)) {
    row_labels <- c(escape_latex(summary_df$algorithm[i]), escape_latex(summary_df$ic_point[i]))
    data_cells <- summary_df[i, -c(1, 2)]
    formatted_cells <- sapply(data_cells, function(cell_content) {
      if (is.na(cell_content)) return("")
      else {
        content_with_break <- str_replace(cell_content, "\\\\n", "\\\\\\\\")
        return(paste0("\\makecell{", content_with_break, "}"))
      }
    })
    full_row <- paste(c(row_labels, formatted_cells), collapse = " & ")
    body_rows <- c(body_rows, paste0(full_row, " \\\\"))
    if (i < nrow(summary_df) && summary_df$algorithm[i] != summary_df$algorithm[i+1]) {
      body_rows <- c(body_rows, "\\addlinespace")
    }
  }

  footer <- c("\\bottomrule", "\\end{tabular}", "\\end{sideways}", "\\end{table*}")
  full_latex_code <- paste(c(header, body_rows, footer), collapse = "\n")
  return(full_latex_code)
}

## fe ----
create_bias_summary_fe <- function(data, rounding_digits = 2) {
  ic_point_original_levels <- levels(data$ic_point)
  ic_mechanism_original_levels <- levels(data$ic_mechanism)

  processed_data <- data %>%
    mutate(
      algorithm = case_when(algorithm == "pam" ~ "PAM", TRUE ~ str_to_title(algorithm)),
      ic_point = str_replace(str_to_title(ic_point), "Adjustment", "Adj."),
      ic_mechanism = str_to_title(ic_mechanism),
      problem_label = case_when(
        problem == "sim_pexp" ~ "Piecewise Exponential DGP",
        problem == "sim_weibull" ~ "Weibull DGP",
        TRUE ~ problem
      )
    )

  wide_table <- processed_data %>%
    pivot_wider(
      id_cols = c(algorithm, ic_point),
      names_from = c(problem_label, ic_mechanism),
      values_from = "bias"
    )

  problems <- c("Piecewise Exponential DGP", "Weibull DGP")
  capitalized_mechanisms <- str_to_title(ic_mechanism_original_levels)
  ordered_data_cols <- as.vector(t(outer(problems, capitalized_mechanisms, paste, sep = "_")))
  desired_col_order <- c("algorithm", "ic_point", ordered_data_cols)

  final_table <- wide_table %>%
    select(any_of(desired_col_order)) %>%
    mutate(across(where(is.numeric), ~ round(., rounding_digits))) %>%
    arrange(
      factor(algorithm, levels = c("PAM", "Cox", "Weibull")),
      factor(ic_point, levels = str_replace(str_to_title(ic_point_original_levels), "Adjustment", "Adj."))
    )

  return(final_table)
}

create_coverage_summary_fe <- function(data, rounding_digits = 2) {
  ic_point_original_levels <- levels(data$ic_point)
  ic_mechanism_original_levels <- levels(data$ic_mechanism)

  sprintf_format <- paste0("%.", rounding_digits, "f\\n(%.", rounding_digits, "f; %.", rounding_digits, "f)")

  processed_data <- data %>%
    mutate(
      algorithm = case_when(algorithm == "pam" ~ "PAM", TRUE ~ str_to_title(algorithm)),
      ic_point = str_replace(str_to_title(ic_point), "Adjustment", "Adj."),
      ic_mechanism = str_to_title(ic_mechanism),
      problem_label = case_when(
        problem == "sim_pexp" ~ "Piecewise Exponential DGP",
        problem == "sim_weibull" ~ "Weibull DGP",
        TRUE ~ problem
      ),
      coverage_formatted = sprintf(sprintf_format, coverage, coverage_lower, coverage_upper)
    )

  wide_table <- processed_data %>%
    pivot_wider(
      id_cols = c(algorithm, ic_point),
      names_from = c(problem_label, ic_mechanism),
      values_from = "coverage_formatted"
    )

  problems <- c("Piecewise Exponential DGP", "Weibull DGP")
  capitalized_mechanisms <- str_to_title(ic_mechanism_original_levels)
  ordered_data_cols <- as.vector(t(outer(problems, capitalized_mechanisms, paste, sep = "_")))
  desired_col_order <- c("algorithm", "ic_point", ordered_data_cols)

  final_table <- wide_table %>%
    select(any_of(desired_col_order)) %>%
    arrange(
      factor(algorithm, levels = c("PAM", "Cox", "Weibull")),
      factor(ic_point, levels = str_replace(str_to_title(ic_point_original_levels), "Adjustment", "Adj."))
    )

  return(final_table)
}


generate_bias_latex_fe <- function(summary_df, rounding_digits = 2) {
  escape_latex <- function(text) {
    if (is.na(text)) return("")
    text %>%
      str_replace_all("&", "\\\\&") %>% str_replace_all("%", "\\\\%") %>%
      str_replace_all("_", "\\\\_")
  }

  sprintf_format <- paste0("%.", rounding_digits, "f")

  data_cols <- names(summary_df)[-c(1, 2)]
  mechanism_headers <- unique(sapply(strsplit(data_cols, "_"), tail, 1))
  header_row <- paste(
    "Algorithm & \\shortstack[l]{Estimation \\\\ Time Point}",
    paste(mechanism_headers, collapse = " & "),
    paste(mechanism_headers, collapse = " & "),
    sep = " & "
  )

  header <- c(
    "\\begin{table*}[htbp]", "\\centering", "\\caption{\\captionicfebias}",
    "\\label{tab:sim-ic-fe-bias}", "\\begin{sideways}", "\\begin{tabular}{llcccccc}",
    "\\toprule", "& & \\multicolumn{3}{c}{Piecewise Exponential DGP} & \\multicolumn{3}{c}{Weibull DGP} \\\\",
    "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
    paste0(header_row, " \\\\"),
    "\\midrule"
  )

  body_rows <- c()
  for (i in 1:nrow(summary_df)) {
    row_labels <- c(escape_latex(summary_df$algorithm[i]), escape_latex(summary_df$ic_point[i]))
    numeric_values <- summary_df[i, -c(1, 2)]
    formatted_values <- sapply(numeric_values, function(val) if (is.na(val)) "" else sprintf(sprintf_format, val))
    full_row <- paste(c(row_labels, formatted_values), collapse = " & ")
    body_rows <- c(body_rows, paste0(full_row, " \\\\"))
    if (i < nrow(summary_df) && summary_df$algorithm[i] != summary_df$algorithm[i + 1]) {
      body_rows <- c(body_rows, "\\addlinespace")
    }
  }

  footer <- c("\\bottomrule", "\\end{tabular}", "\\end{sideways}", "\\end{table*}")
  full_latex_code <- paste(c(header, body_rows, footer), collapse = "\n")
  return(full_latex_code)
}

generate_coverage_latex_fe <- function(summary_df) {
  escape_latex <- function(text) {
    if (is.na(text)) return("")
    text %>%
      str_replace_all("&", "\\\\&") %>% str_replace_all("%", "\\\\%") %>%
      str_replace_all("_", "\\\\_")
  }

  data_cols <- names(summary_df)[-c(1, 2)]
  mechanism_headers <- unique(sapply(strsplit(data_cols, "_"), tail, 1))
  header_row <- paste(
    "Algorithm & \\shortstack[l]{Estimation \\\\ Time Point}",
    paste(mechanism_headers, collapse = " & "),
    paste(mechanism_headers, collapse = " & "),
    sep = " & "
  )

  header <- c(
    "\\begin{table*}[htbp]", "\\centering", "\\caption{\\captionicfecoverage}",
    "\\label{tab:sim-ic-fe-coverage}", "\\begin{sideways}",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\begin{tabular}{llcccccc}",
    "\\toprule", "& & \\multicolumn{3}{c}{Piecewise Exponential DGP} & \\multicolumn{3}{c}{Weibull DGP} \\\\",
    "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
    paste0(header_row, " \\\\"),
    "\\midrule"
  )

  body_rows <- c()
  for (i in 1:nrow(summary_df)) {
    row_labels <- c(escape_latex(summary_df$algorithm[i]), escape_latex(summary_df$ic_point[i]))
    data_cells <- summary_df[i, -c(1, 2)]
    formatted_cells <- sapply(data_cells, function(cell_content) {
      if (is.na(cell_content)) return("")
      else {
        content_with_break <- str_replace(cell_content, "\\\\n", "\\\\\\\\")
        return(paste0("\\makecell{", content_with_break, "}"))
      }
    })
    full_row <- paste(c(row_labels, formatted_cells), collapse = " & ")
    body_rows <- c(body_rows, paste0(full_row, " \\\\"))
    if (i < nrow(summary_df) && summary_df$algorithm[i] != summary_df$algorithm[i + 1]) {
      body_rows <- c(body_rows, "\\addlinespace")
    }
  }

  footer <- c("\\bottomrule", "\\end{tabular}", "\\end{sideways}", "\\end{table*}")
  full_latex_code <- paste(c(header, body_rows, footer), collapse = "\n")
  return(full_latex_code)
}
