library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(mvtnorm)
library(rlang)
library(survival)
library(icenReg)
library(mgcv)
library(pammtools)

read_file <- function(a_infile, header="auto",sep="auto",fill=FALSE, ...){
	tblInRaw <- fread(a_infile, header=header, sep=sep)
	tblIn <- as.data.frame(tblInRaw,fill=fill)
	return(tblIn)
}


sim_wrapper <- function(
  data,
  job,
  formula = NULL,
  n = 500,
  time_grid = seq(0, 10, by = 0.05),
  ic = TRUE,
  ic_mechanism = c("beta", "uniform"),
  round = NULL) {

  # create data set with covariates
  df <- tibble::tibble(x1 = rbinom(n, 1, 0.8), x2 = runif(n, 0, 6))

  if(is.null(formula)) {
    formula <- as.formula("~ -3.5 + dgamma(t, 8, 2)*6 - 1.3*x1 + sqrt(x2)")
  } else {
     formula <- as.formula(formula)
  }

  ndf <- sim_pexp(formula = formula, data = df, cut = time_grid)

  if(ic) {
    ndf <- add_interval_censoring(ndf, max_time = max(ndf$time), visits_range = c(1, 10), ic_mechanism = ic_mechanism, round = round)
  }

  out <- list(ndf, formula)

  return(out)
}


coverage_wrapper_pam <- function(
  data,
  job,
  instance,
  bs = "ps",
  k = 10,
  # ic_point = c("start", "mid", "end", "oracle")
  ic_point = c("mid", "end", "oracle")) {

  df <- instance[[1]]
  formula <- instance[[2]]

  if(!ic_point %in% c("start", "mid", "end", "oracle"))
    stop("ic_point must be 'mid', 'end' or 'oracle'!")
    # stop("ic_point must be 'start', 'mid', 'end' or 'oracle'!")

  if(ic_point == "mid") {
    df <- df %>%
      mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    # ic_point == "start" ~ "Surv(time_ic_start, status_ic) ~ x1 + x2",
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
    ic_point == "end" ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
    ic_point == "oracle" ~ "Surv(time, status) ~ x1 + x2",
    TRUE ~ "NA") %>%
    as.formula()

  ped <- as_ped(
    data    = df,
    formula = formula_ped,
    id      = "id")

  formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + x1 + s(x2)")

  mod <- bam(
    formula = as.formula(formula_mod),
    data    = ped,
    family  = poisson(),
    offset  = offset,
    method  = "fREML",
    discrete = TRUE)

  # create new data set
  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- make_newdata(
    ped,
    tend = sort(unique(ped$tend)),
    x1 = median(df$x1),
    x2 = mean(df$x2)) %>%
    mutate(
      true_hazard = exp(eval(parse(text = formula), envir = pick(everything()))),
      true_cumu = cumsum(intlen * true_hazard),
      true_surv = exp(-true_cumu)) %>%
    add_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    add_cumu_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    add_surv_prob(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
    mutate(
      hazard = (true_hazard >= ci_lower) & (true_hazard <= ci_upper),
      cumu = (true_cumu >= cumu_lower) & (true_cumu <= cumu_upper),
      surv =  (true_surv >= surv_lower) & (true_surv <= surv_upper)) %>%
    select(hazard, cumu, surv) %>%
    summarize_all(mean)

    return(nd)
}


coverage_wrapper_cox <- function(
  data,
  job,
  instance,
  formula = NULL,
  ic_point = c("mid", "end", "oracle")) {

  if(!ic_point %in% c("start", "mid", "end", "oracle"))
    stop("ic_point must be 'mid', 'end' or 'oracle'!")

  df <- instance[[1]]
  formula <- instance[[2]]

  # fit regular coxph model
  formula_cox <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + pspline(x2)",
    ic_point == "end" ~ "Surv(time_ic_stop, status_ic) ~ x1 + pspline(x2)",
    ic_point == "oracle" ~ "Surv(time, status) ~ x1 + pspline(x2)",
    TRUE ~ "NA") %>%
    as.formula()

  mod <- coxph(formula = formula_cox, data = df)

  # create new dataset exactly as specified
  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
    ic_point == "end" ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
    ic_point == "oracle" ~ "Surv(time, status) ~ x1 + x2",
    TRUE ~ "NA") %>%
    as.formula()

  ped <- as_ped(
    data    = df,
    formula = formula_ped,
    id      = "id")

  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- make_newdata(
    ped,
    tend = sort(unique(ped$tend)),
    x1 = median(df$x1),
    x2 = mean(df$x2)
  ) %>%
  mutate(
    true_hazard = exp(eval(parse(text = formula), envir = pick(everything()))),
    true_cumu = cumsum(intlen * true_hazard),
    true_surv = exp(-true_cumu)
  )

  # get baseline estimates and CI via survfit
  newdata_cox <- data.frame(x1=median(df$x1), x2=mean(df$x2))
  sf <- survfit(mod, newdata = newdata_cox)

  # Extract baseline values at matching times
  est_surv <- approx(sf$time, sf$surv, nd$tend, rule=2)$y
  surv_upper <- approx(sf$time, sf$upper, nd$tend, rule=2)$y
  surv_lower <- approx(sf$time, sf$lower, nd$tend, rule=2)$y

  est_cumu <- -log(est_surv)
  cumu_upper <- -log(surv_lower)
  cumu_lower <- -log(surv_upper)

  baseline_hazard <- c(est_cumu[1]/nd$tend[1], diff(est_cumu)/diff(nd$tend))

  # Linear predictor and its SE
  lp_pred <- predict(mod, newdata = newdata_cox, type="lp", se.fit=TRUE)
  lp <- lp_pred$fit
  lp_se <- lp_pred$se.fit

  # hazard, cumu hazard, and survival adjusted by linear predictor
  nd <- nd %>%
    mutate(
      est_hazard = baseline_hazard * exp(lp), # BUG: no uncertainty in baseline_hazard!!!
      est_cumu = est_cumu,
      est_surv = est_surv,
      # approximate hazard CI by propagating lp_se (CI on LP)
      hazard_upper = baseline_hazard * exp(lp + 1.96*lp_se),
      hazard_lower = baseline_hazard * exp(lp - 1.96*lp_se),
      cumu_upper = cumu_upper,
      cumu_lower = cumu_lower,
      surv_upper = surv_upper,
      surv_lower = surv_lower
    ) %>%
    mutate(
      hazard_cov = (true_hazard >= hazard_lower) & (true_hazard <= hazard_upper),
      cumu_cov = (true_cumu >= cumu_lower) & (true_cumu <= cumu_upper),
      surv_cov = (true_surv >= surv_lower) & (true_surv <= surv_upper)
    ) %>%
    select(hazard_cov, cumu_cov, surv_cov) %>%
    summarize_all(mean)

  return(nd)
}


add_interval_censoring <- function(data, max_time, visits_range = c(1, 10), ic_mechanism = c("beta", "uniform"), round = NULL) {

  # check necessary columns
  required_cols <- c("id", "time", "status")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: id, time, and status.")
  }
  if (!ic_mechanism %in% c("beta", "uniform")) {
    stop("ic_mechanism must be 'beta' or 'uniform'.")
  }

  interval_censor_individual_beta <- function(event_time, event_status, visits_range, max_time, round) {
    v <- sample(visits_range[1]:visits_range[2], 1)
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

  interval_censor_individual_uniform <- function(event_time, event_status, visits_range, max_time, round) {

    # tbd!!!

  }

  # apply to each individual
  interval_fun <- if (ic_mechanism == "beta") {
    interval_censor_individual_beta
  } else {
    interval_censor_individual_uniform
  }

  interval_data <- t(mapply(interval_fun,
                            event_time = data$time,
                            event_status = data$status,
                            MoreArgs = list(visits_range = visits_range,
                                            max_time = max_time,
                                            round = round)))

  interval_data_df <- as.data.frame(interval_data)

  # return the original data with interval-censored data columns
  final_data <- dplyr::bind_cols(data, interval_data_df)

  return(final_data)
}

