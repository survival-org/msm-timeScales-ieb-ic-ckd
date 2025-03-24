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
library(knitr)
library(ggplot2)


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


wrapper_pam <- function(
  data,
  job,
  instance,
  bs = "ps",
  k = 10,
  ic_point = c("mid", "end", "mid_end", "true_time", "oracle")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  if(ic_point == "mid") {
    df <- df %>%
      mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
    ic_point %in% c("end", "mid_end") ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
    ic_point %in% c("true_time", "oracle") ~ "Surv(time, status) ~ x1 + x2",
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
  if(ic_point == "oracle") {
    formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + x1 + sqrt(x2)")
  } else {
    formula_mod <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + x1 + s(x2)")
  }

  if(ic_point != "mid_end") {
    mod <- bam(
      formula = as.formula(formula_mod),
      data    = ped,
      family  = poisson(),
      offset = offset,
      method  = "fREML",
      discrete = TRUE)
  } else {
    mod <- bam(
      formula = as.formula(formula_mod),
      data    = ped %>% mutate(tend = tmid) %>% select(-tmid),
      family  = poisson(),
      offset = offset,
      method  = "fREML",
      discrete = TRUE)
  }

  # create new data set
  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- make_newdata(
    ped,
    tend = sort(unique(ped$tend)),
    x1 = median(df$x1),
    x2 = mean(df$x2))

  if(ic_point != "mid_end") {
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
  } else {
    nd <- nd %>%
      mutate(tmid = (tstart + tend) / 2) %>%
      mutate(
        loghazard_true = eval(parse(text = formula), envir = pick(everything())),
        hazard_true = exp(loghazard_true),
        cumu_true = cumsum(intlen * hazard_true),
        surv_true = exp(-cumu_true)) %>%
      rename(tend_original = tend, tend = tmid) %>%
      add_hazard(mod, se_mult = qnorm(0.975), ci_type = "default", type = "link") %>%
      rename(loghazard = hazard, loghazard_lower = ci_lower, loghazard_upper = ci_upper) %>%
      add_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
      rename(hazard_lower = ci_lower, hazard_upper = ci_upper) %>%
      add_cumu_hazard(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
      rename(cumu = cumu_hazard) %>%
      add_surv_prob(mod, se_mult = qnorm(0.975), ci_type = "default") %>%
      rename(surv = surv_prob) %>%
      rename(tmid = tend, tend = tend_original) %>%
      mutate(
        loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
        hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
        cumu_cov = as.integer((cumu_true >= cumu_lower) & (cumu_true <= cumu_upper)),
        surv_cov =  as.integer((surv_true >= surv_lower) & (surv_true <= surv_upper)))
  }

  nd <- nd %>%
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov)

  return(nd)
}


wrapper_cox <- function(
  data,
  job,
  instance,
  formula = NULL,
  ic_point = c("mid", "end", "true_time", "oracle")) {

  ic_point <- match.arg(ic_point)

  df <- instance[[1]]
  formula <- instance[[2]]

  if(ic_point == "mid") {
    df <- df %>%
      mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
    ic_point == "end" ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
    ic_point %in% c("true_time", "oracle") ~ "Surv(time, status) ~ x1 + x2",
    TRUE ~ NA_character_) %>%
    as.formula()

  ped <- as_ped(
    data    = df,
    formula = formula_ped,
    id      = "id")

  # fit regular coxph model
  formula_mod <- update(formula_ped, . ~ x1 + pspline(x2))

  mod <- coxph(formula = formula_mod, data = df)

  # create new dataset exactly as specified
  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- make_newdata(
    ped,
    tend = sort(unique(ped$tend)),
    x1 = median(df$x1),
    x2 = mean(df$x2)
  ) %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)
    )

  # get baseline estimates and CI via survfit
  newdata_cox <- data.frame(x1=median(df$x1), x2=mean(df$x2))
  sf <- survfit(mod, newdata = newdata_cox)

  # Extract baseline values at matching times
  surv <- approx(sf$time, sf$surv, nd$tend, rule=2)$y
  surv_upper <- approx(sf$time, sf$upper, nd$tend, rule=2)$y
  surv_lower <- approx(sf$time, sf$lower, nd$tend, rule=2)$y

  cumu <- -log(surv)
  cumu_upper <- -log(surv_lower)
  cumu_lower <- -log(surv_upper)

  baseline_hazard <- c(cumu[1]/nd$tend[1], diff(cumu)/diff(nd$tend))

  # Linear predictor and its SE
  lp_pred <- predict(mod, newdata = newdata_cox, type="lp", se.fit=TRUE)
  lp <- lp_pred$fit
  lp_se <- lp_pred$se.fit

  # hazard, cumu hazard, and survival adjusted by linear predictor
  nd <- nd %>%
    mutate(
      loghazard = NA_real_,
      hazard = baseline_hazard * exp(lp), # BUG: no uncertainty in baseline_hazard!!!
      cumu = cumu,
      surv = surv,
      # approximate hazard CI by propagating lp_se (CI on LP)
      loghazard_lower = NA_real_,
      loghazard_upper = NA_real_,
      hazard_upper = baseline_hazard * exp(lp + 1.96*lp_se),
      hazard_lower = baseline_hazard * exp(lp - 1.96*lp_se),
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
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov)

  return(nd)
}


wrapper_generalizedGamma <- function(
  data,
  job,
  instance,
  ic_point = c("mid", "end", "true_time", "oracle", "adjustment")) {

  ic_point <- match.arg(ic_point)
saveRDS(instance, "instance.rds")
  df <- instance[[1]]
  formula <- instance[[2]]

  if(ic_point == "mid") {
    df <- df %>% mutate(time_ic_mid = (time_ic_start + time_ic_stop) / 2)
  }

  formula_ped <- case_when(
    ic_point == "mid" ~ "Surv(time_ic_mid, status_ic) ~ x1 + x2",
    ic_point %in% c("end", "adjustment") ~ "Surv(time_ic_stop, status_ic) ~ x1 + x2",
    ic_point %in% c("true_time", "oracle") ~ "Surv(time, status) ~ x1 + x2",
    TRUE ~ NA_character_
  ) %>% as.formula()

  ped <- as_ped(
    data = df,
    formula = formula_ped,
    id = "id")

  if(ic_point != "adjustment") {
    # regular formula
    formula_mod <- update(formula_ped, . ~ x1 + pspline(x2))
  }  else {
    # interval-censored formula
    formula_mod <- Surv(time = time_ic_start, time2 = time_ic_stop, type = "interval2") ~ x1 + pspline(x2)
  }

  # Model fitting
  mod <- tryCatch(
    withCallingHandlers(
      {
        flexsurvreg(formula_mod, data = df, dist = "gengamma",
                    inits = c(mu = log(median(df$time)), sigma = 0.5, Q = 1),
                    method = "BFGS")
      },
      warning = function(w) {
        message("Warning detected in BFGS: Falling back to Nelder-Mead")
        invokeRestart("muffleWarning")  # Prevents warning from propagating further
      }
    ),
    error = function(e) {
      message("Error in BFGS: Falling back to Nelder-Mead")
      flexsurvreg(formula_mod, data = df, dist = "gengamma",
                  inits = c(mu = log(median(df$time)), sigma = 0.5, Q = 1),
                  method = "Nelder-Mead")
    }
  )

  # new data preparation
  formula <- paste(gsub("\\bt\\b", "tend", deparse(formula[[2]])))

  nd <- make_newdata(
    ped,
    tend = sort(unique(ped$tend)),
    x1 = median(df$x1),
    x2 = mean(df$x2)) %>%
    mutate(
      loghazard_true = eval(parse(text = formula), envir = pick(everything())),
      hazard_true = exp(loghazard_true),
      cumu_true = cumsum(intlen * hazard_true),
      surv_true = exp(-cumu_true)
    )

  # predictions from parametric model
  pred <- summary(mod, newdata = data.frame(x1 = median(df$x1), x2 = mean(df$x2)),
                  type = "hazard", t = nd$tend)

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
    select(
      tstart, tend,
      loghazard, hazard, cumu, surv,
      loghazard_true, hazard_true, cumu_true, surv_true,
      loghazard_cov, hazard_cov, cumu_cov, surv_cov)

  return(nd)
}


add_interval_censoring <- function(
  data,
  max_time,
  visits_range = c(1, 10),
  ic_mechanism = c("beta", "uniform"),
  round = NULL) {

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

  interval_censor_individual_uniform <- function(event_time, event_status, visits_range, max_time, round_digits) {
    v <- sample(visits_range[1]:visits_range[2], 1)

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

  # apply to each individual
  interval_fun <- if (ic_mechanism == "beta") {
    interval_censor_individual_beta
  } else {
    interval_censor_individual_uniform
  }

  interval_list <- mapply(interval_fun,
                              event_time = data$time,
                              event_status = data$status,
                              MoreArgs = list(visits_range = visits_range,
                                              max_time = max_time,
                                              round = round),
                              SIMPLIFY = FALSE)

  interval_data <- do.call(rbind, interval_list)

  interval_data_df <- as.data.frame(interval_data)

  # return the original data with interval-censored data columns
  final_data <- bind_cols(data, interval_data_df)

  return(final_data)
}


calc_coverage <- function(data, grouping_vars = NULL, rounding = 3) {

  required_cols <- c(
    grouping_vars,
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

  return(kable(results))
}


calc_rmse <- function(data, grouping_vars = NULL, rounding = 3) {

  required_cols <- c(
    grouping_vars,
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

  return(kable(results))
}


create_linePlot <- function(data, grouping_vars = NULL) {

  required_cols <- c(
    grouping_vars, "tend", "loghazard", "loghazard_true", "job.id"
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  plot_data <- data %>%
    filter(!is.na(loghazard))

  if (!is.null(grouping_vars)) {
    plot_data <- plot_data %>%
      unite("grouping", all_of(grouping_vars), remove = FALSE)
  } else {
    plot_data$grouping <- "All"
  }

  unique_groups <- unique(plot_data$grouping)

  plots <- lapply(unique_groups, function(grp) {
    df_grp <- plot_data %>% filter(grouping == grp)

    ggplot(df_grp, aes(x = tend, y = loghazard)) +
      geom_step(aes(group = job.id), alpha = 0.3) +
      geom_line(aes(y = loghazard_true, col = "truth"), lwd = 1.5) +
      geom_smooth(aes(col = "average estimate"), method = "gam", formula = y ~ s(x), se = FALSE) +
      scale_color_brewer("", palette = "Dark2") +
      xlab("time") +
      ylab("loghazard") +
      ggtitle(paste("Group:", grp)) +
      theme_bw()
  })

  names(plots) <- unique_groups

  return(plots)
}
