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

# data prep ----
wrapper_sim <- function(
  data,
  job,
  formulas_list,
  terminal_states,
  cut,
  n = 2500,
  round = 2,
  cens_type = c("none", "right", "interval", "left"),
  cens_dist = c("weibull", "exponential", "lognormal", "uniform"),
  cens_params = NULL){

  data <- data.frame(
    id          = seq_len(n),
    x1          = rbinom(n, 1, 0.5),
    from        = 0L,
    t           = 0 # “calendar” clock
  )

  if(is.list(terminal_states)) {
    terminal_states <- unlist(terminal_states)
  }

  if(!is.null(cut)) {
    cut <- unlist(cut)
  }

  df <- sim_pexp_msm(formulas_list, data, cut, terminal_states, round = round, add_counterfactuals = FALSE)

  if(cens_type != "none") {

    if (is.list(cens_params)) {
      cens_params <- unlist(cens_params)
    }

    df <- df %>%
      add_censoring(type = cens_type, distribution = cens_dist, parameters = unlist(cens_params), round = round)
  }

  ped <- as_ped_multistate(
    data       = df %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

  ped <- ped %>%
    add_timeScales() %>%
    add_transVars()

  out <- list(
    formulas_list = formulas_list,
    df = df,
    ped = ped
  )

  return(out)
}


wrapper_baselineHazards <- function(
  data,
  job,
  instance,
  formula) {

  formulas_list <- instance$formulas_list
  ped <- instance$ped

  # TBD
}


wrapper_fixedEffects <- function(
  data,
  job,
  instance,
  formula) {

  formulas_list <- instance$formulas_list
  ped <- instance$ped

  # Extract transitions from the formulas_list
  transitions <- map_chr(formulas_list, ~ {
    from <- .x$from
    to   <- .x$to
    paste0(from, "->", to)
  })

  # Extract all beta_1 coefficients from formulas_list
  beta_1_vars <- formulas_list %>%
    map(~ {
      formula <- .x$formula
      terms <- all.vars(formula)
      beta_1_term <- terms[grepl("^beta_1_", terms)]
      beta_1_term
    }) %>%
    unlist()

  # Evaluate the extracted coefficients
  beta_1_true <- beta_1_vars %>%
    map(~ eval(parse(text = .x), envir = ped)) %>%
    set_names(beta_1_vars)

  mod <- bam(as.formula(formula)
                      , data = ped
                      , family=poisson()
                      , offset=offset
                      , discrete = T
                      , method = "fREML"
  )
  summary_mod <- summary(mod)

  # Adjust coefficients and standard errors
  beta_est <- map_dbl(beta_1_vars, ~ {
    if (.x == beta_1_vars[1]) {
      summary_mod$p.coeff["x1"]
    } else {
      transition_name <- transitions[which(beta_1_vars == .x)]
      summary_mod$p.coeff[paste0("transition", transition_name, ":", "x1")]
    }
  })

  beta_se <- map_dbl(beta_1_vars, ~ {
    if (.x == beta_1_vars[1]) {
      summary_mod$se["x1"]
    } else {
      transition_name <- transitions[which(beta_1_vars == .x)]
      summary_mod$se[paste0("transition", transition_name, ":", "x1")]
    }
  })

  # Calculate confidence intervals
  ci_lower <- beta_est - 1.96 * beta_se
  ci_upper <- beta_est + 1.96 * beta_se

  # Determine coverage
  beta_cov <- map2_dbl(beta_1_true, seq_along(beta_1_vars), ~ {
    if (.y == 1) {
      ifelse(.x >= ci_lower[.y] && .x <= ci_upper[.y], 1, 0)
    } else {
      transition_name <- transitions[.y]
      ifelse(.x >= ci_lower[.y] && .x <= ci_upper[.y], 1, 0)
    }
  })

  coef_df <- tibble(
    variable   = beta_1_vars,                # "beta_1_01", "beta_1_04", ...
    true_value = unlist(beta_1_true),        # the “ground‐truth” from your sim
    estimate   = beta_est,                    # fitted estimate
    std_error  = beta_se,                     # its SE
    ci_lower   = ci_lower,                    # 95% CI lower
    ci_upper   = ci_upper,                    # 95% CI upper
    coverage   = beta_cov                     # 1 if true_value ∈ [ci_lower, ci_upper]
  )

  coef_df <- coef_df %>%
    mutate(
      estimate  = round(estimate,  4),
      std_error = round(std_error, 4),
      ci_lower  = round(ci_lower,  4),
      ci_upper  = round(ci_upper,  4)
    )

  return(coef_df)
}


#' Draw random numbers from piece-wise exponential distribution.
#'
#' This is a copy of the same function from \code{rpexp} from package
#' \pkg{msm}.
#' Copied here to reduce dependencies.
#'
#' @inheritParams msm::rpexp
#' @importFrom stats rexp
#'
#' @keywords internal
rpexp_msm <- function (n = 1, rate = 1, t = 0) {

    if (length(t) != length(rate))
        stop("length of t must be equal to length of rate")
    # if (!isTRUE(all.equal(0, t[1])))
    #     stop("first element of t should be 0")
    if (is.unsorted(t))
        stop("t should be in increasing order")
    if (length(n) > 1)
        n <- length(n)
    if (n == 0)
        return(numeric(0))
    if (length(rate) == 1)
        return(rexp(n, rate))
    ret <- numeric(n)
    left <- 1:n
    for (i in seq_along(rate)) {
        re <- rexp(length(left), rate[i])
        r <- t[i] + re
        success <- if (i == length(rate))
            seq_along(left)
        else which(r < t[i + 1])
        ret[left[success]] <- r[success]
        left <- setdiff(left, left[success])
        if (length(left) == 0)
            break
    }
    ret
}


#' Simulate survival times from the piece-wise exponential distribution
#'
#' @param formula An extended formula that specifies the linear predictor.
#' If you want to include a smooth baseline
#' or time-varying effects, use \code{t} within your formula as
#' if it was a covariate in the data, although it is not and should not
#' be included in the \code{data} provided to \code{sim_pexp}. See examples
#' below.
#'
#' @param data A data set with variables specified in \code{formula}.
#' @param cut A sequence of time-points starting with 0.
#' @import dplyr
#' @import Formula
#' @importFrom lazyeval f_eval
#' @importFrom tidyr replace_na
#' @examples
#' library(survival)
#' library(dplyr)
#' library(pammtools)
#'
#' # set number of observations/subjects
#' n <- 250
#' # create data set with variables which will affect the hazard rate.
#' df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
#'  as_tibble()
#' # the formula which specifies how covariates affet the hazard rate
#' f0 <- function(t) {
#'  dgamma(t, 8, 2) *6
#' }
#' form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)
#' set.seed(24032018)
#' sim_df <- sim_pexp(form, df, 1:10)
#' head(sim_df)
#' plot(survfit(Surv(time, status)~1, data = sim_df ))
#'
#' # for control, estimate with Cox PH
#' mod <- coxph(Surv(time, status) ~ x1 + pspline(x2), data=sim_df)
#' coef(mod)[1]
#' layout(matrix(1:2, nrow=1))
#' termplot(mod, se = TRUE)
#'
#' # and using PAMs
#' layout(1)
#' ped <- sim_df %>% as_ped(Surv(time, status)~., max_time=10)
#' library(mgcv)
#' pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
#' coef(pam)[2]
#' plot(pam, page=1)
#'
#'\dontrun{
#' # Example 2: Functional covariates/cumulative coefficients
#' # function to generate one exposure profile, tz is a vector of time points
#' # at which TDC z was observed
#' rng_z = function(nz) {
#'   as.numeric(arima.sim(n = nz, list(ar = c(.8, -.6))))
#' }
#' # two different exposure times  for two different exposures
#' tz1 <- 1:10
#' tz2 <- -5:5
#' # generate exposures and add to data set
#' df <- df %>%
#'   add_tdc(tz1, rng_z) %>%
#'   add_tdc(tz2, rng_z)
#' df
#'
#' # define tri-variate function of time, exposure time and exposure z
#' ft <- function(t, tmax) {
#'   -1*cos(t/tmax*pi)
#' }
#' fdnorm <- function(x) (dnorm(x,1.5,2)+1.5*dnorm(x,7.5,1))
#' wpeak2 <- function(lag) 15*dnorm(lag,8,10)
#' wdnorm <- function(lag) 5*(dnorm(lag,4,6)+dnorm(lag,25,4))
#' f_xyz1 <- function(t, tz, z) {
#'   ft(t, tmax=10) * 0.8*fdnorm(z)* wpeak2(t - tz)
#' }
#' f_xyz2 <- function(t, tz, z) {
#'   wdnorm(t-tz) * z
#' }
#'
#' # define lag-lead window function
#' ll_fun <- function(t, tz) {t >= tz}
#' ll_fun2 <- function(t, tz) {t - 2 >= tz}
#' # simulate data with cumulative effect
#' sim_df <- sim_pexp(
#'   formula = ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)|
#'      fcumu(t, tz1, z.tz1, f_xyz=f_xyz1, ll_fun=ll_fun) +
#'      fcumu(t, tz2, z.tz2, f_xyz=f_xyz2, ll_fun=ll_fun2),
#'   data = df,
#'   cut = 0:10)
#'}
#' @export
sim_pexp <- function(formula, data, cut) {

  data <- data %>%
    mutate(
      id     = row_number(),
      time   = max(cut),
      status = 1)

  # extract formulas for different components
  Form <- Formula(formula)
  f1   <- formula(Form, rhs = 1)
  # later more sophisticated checks + could be used to map over all rhs
  # formulae, check what type of evaluation is needed and return ETAs for
  # each part of the formula separated by |, such that model estimation may
  # be checked for individuals terms/parts
  if (length(Form)[2] > 1) {
    f2  <- formula(Form, rhs = 2)
  } else {
    f2 <- NULL
  }

  # construct eta for time-constant part
  ped  <- split_data(
      formula = Surv(time, status)~.,
      data    = select_if(data, is_atomic),
      cut     = cut,
      id      = "id") %>%
    rename("t" = "tstart") %>%
    mutate(rate = exp(f_eval(f1, .)))

  # construct eta for time-dependent part
  if (!is.null(f2)) {
    terms_f2  <- terms(f2, specials = "fcumu")
    f2_ev     <- list()
    f2_tl <- attr(terms_f2, "term.labels")
    for (i in seq_along(f2_tl)) {
      f2_ev[[i]] <- eval(expr = parse(text = f2_tl[[i]]), envir = .GlobalEnv)
    }
    ll_funs   <- map(f2_ev, ~.x[["ll_fun"]])
    tz_vars   <- map_chr(f2_ev, ~.x[["vars"]][1])
    cumu_funs <- map(f2_ev, ~.x[["f_xyz"]])
    names(tz_vars) <- names(ll_funs) <- names(cumu_funs) <- tz_vars
    z_form <- list("eta_", map_chr(f2_ev, ~.x[["vars"]][2])) %>%
      reduce(paste0, collapse = "+") %>% paste0("~", .) %>% as.formula()

    df2 <- map(f2_ev, function(fc) eta_cumu(data = data, fc, cut = cut))
    suppressMessages(
      ped <- ped %>%
        left_join(reduce(df2, full_join))
    )
    ped <- ped %>%
      mutate_at(vars(contains("eta_")), replace_na, 0) %>%
      group_by(.data$id, .data$t) %>%
      mutate(eta_z = !!rlang::get_expr(z_form)) %>%
      mutate(rate = .data$rate * exp(.data$eta_z))
  } else {
    tz_vars <- NULL
  }

  sim_df <- ped %>%
    group_by(id) %>%
    summarize(time = pammtools:::rpexp(rate = .data$rate, t = .data$t)) %>%
    mutate(
      status = 1L * (.data$time <= max(cut)),
      time   = pmin(.data$time, max(cut)))

  suppressMessages(
    sim_df <- sim_df %>%
      left_join(select(data, -all_of(c("time", "status"))))
  )

  attr(sim_df, "id_var")     <- "id"
  attr(sim_df, "time_var")   <- "time"
  attr(sim_df, "status_var") <- "status"
  attr(sim_df, "tz_var")     <- tz_vars
  attr(sim_df, "cens_value") <- 0
  attr(sim_df, "breaks")     <- cut
  attr(sim_df, "tz")         <- imap(tz_vars, ~select(sim_df, all_of(.x)) %>%
    pull(.x) %>% unique()) %>% flatten()
  if (exists("ll_funs")) attr(sim_df, "ll_funs") <- ll_funs
  if (exists("cumu_funs")) attr(sim_df, "cumu_funs") <- cumu_funs
  attr(sim_df, "id_n") <- sim_df %>% pull("time") %>%
    pmin(max(cut)) %>%
    map_int(findInterval, vec = cut, left.open = TRUE, rightmost.closed = TRUE)
  attr(sim_df, "id_tseq") <- attr(sim_df, "id_n") %>%
    map(seq_len) %>% unlist()
  attr(sim_df, "id_tz_seq") <- rep(seq_along(pull(sim_df, id)),
    times = attr(sim_df, "id_n"))
  attr(sim_df, "sim_formula") <- formula

  class(sim_df) <- c("sim_df", class(unped(sim_df)))

  if (any(!map_lgl(sim_df, is_atomic))) {
    class(sim_df) <- c("nested_fdf", class(sim_df))
  }

  sim_df

}


#' Add time-dependent covariate to a data set
#'
#' Given a data set in standard format (with one row per subject/observation),
#' this function adds a column with the specified exposure time points
#' and a column with respective exposures, created from \code{rng_fun}.
#' This function should usually only be used to create data sets passed
#' to \code{\link[pammtools]{sim_pexp}}.
#'
#' @inheritParams sim_pexp
#' @param tz A numeric vector of exposure times (relative to the
#' beginning of the follow-up time \code{t})
#' @param rng_fun A random number generating function that creates
#' the time-dependent covariates at time points \code{tz}.
#' First argument of the function should be \code{n}, the number of
#' random numbers to generate. Within \code{add_tdc}, \code{n} will be set
#' to \code{length(tz)}.
#' @param ... Currently not used.
#' @import dplyr
#' @importFrom rlang eval_tidy :=
#' @importFrom purrr map
#' @export
add_tdc <- function(data, tz, rng_fun, ...) {

  tz      <- enquo(tz)
  nz      <- length(eval_tidy(tz))
  name_tz <- quo_name(tz)
  z_var   <- paste0("z.", name_tz)

  data %>%
    mutate(
      !!name_tz := map(seq_len(n()), ~ !!tz),
      !!z_var   := map(seq_len(n()), ~ rng_fun(nz = nz))) %>%
    as_tibble()

}


#' A formula special used to handle cumulative effect specifications
#'
#' Can be used in the second part of the formula specification provided
#' to \code{\link[pammtools]{sim_pexp}} and should only be used in this
#' context.
#'
#' @importFrom purrr map
#' @export
#' @keywords internal
fcumu <- function(..., by = NULL, f_xyz, ll_fun) {

  vars   <- as.list(substitute(list(...)))[-1] %>%
    map(~as.character(.x)) %>%
    unlist()
  vars <- vars[vars != "t"]

  list(
    vars   = vars,
    f_xyz  = f_xyz,
    ll_fun = ll_fun)

}


#' @import dplyr
#' @importFrom tidyr unnest
#' @importFrom rlang sym :=
#' @keywords internal
eta_cumu <- function(data, fcumu, cut, ...) {

  vars   <- fcumu$vars
  f_xyz  <- fcumu$f_xyz
  ll_fun <- fcumu$ll_fun
  eta_name <- paste0("eta_", vars[2])
  comb_df <- combine_df(
    data.frame(t = cut),
    select(data, one_of("id", vars)))
  comb_df <- comb_df %>% unnest(cols = -one_of("id"))
  comb_df %>%
    group_by(.data$id, .data$t) %>%
    mutate(
      LL = ll_fun(t, !!sym(vars[1])) * 1,
      delta = c(mean(abs(diff(!!sym(vars[1])))), abs(diff(!!sym(vars[1]))))) %>%
    ungroup() %>%
    filter(.data$LL != 0) %>%
    group_by(.data$id, .data$t) %>%
    summarize(!!eta_name :=
      sum(.data$delta * f_xyz(.data$t, .data[[vars[1]]], .data[[vars[2]]])))

}


#' Simulate data for competing risks scenario
#'
#'
#' @keywords internal
sim_pexp_cr_msm <- function(formula, data, from, to, cut, round = NULL) {

  # Validate 'from' and 'to'
  if (length(from) != 1) stop("Argument 'from' must be a single state.")

  # Formula extends the base class formula by allowing for multiple responses and multiple parts of regressors
  Form    <- Formula(formula)
  # Extract the right handside of the Formula
  F_rhs   <- attr(Form, "rhs")
  l_rhs   <- length(F_rhs)
  if (length(to) != l_rhs) stop(sprintf(
    "Argument 'to' must have length %d to match the %d hazards in formula.",
    l_rhs, l_rhs
  ))
  seq_rhs <- seq_len(l_rhs)

  if (!("id" %in% names(data))) {
    data$id <- 1:(nrow(data))
  }

  if (!("t" %in% names(data))) {
    data$t <- 0
  }

  data <- data %>%
    mutate(
      time_start = t,
      time_end   = max(cut),
      status = 1
    )

  # Update clocks: t_until_<from> always = t
  data[[paste0("t_until_", from)]] <- data$t

  # construct eta for time-constant part
  # offset (the log of the duration during which the subject was under risk in that interval)
  ped  <- split_data(
    formula = Surv(time_start, time_end, status) ~ .,
    data    = select_if(data, is_atomic),
    cut     = cut,
    id      = "id")

  # find all t_until_ columns
  t_until_cols <- grep("^t_until_", names(ped), value = TRUE)
  # for each, create matching t_<k> = ped$tend - ped[[t_until_col]]
  for (col in t_until_cols) {
    k      <- sub("^t_until_", "", col)
    outcol <- paste0("t_", k)
    ped[[outcol]] <- with(ped,
      ifelse(
        is.na(ped[[col]]),
        NA_real_,
        ped$tend - ped[[col]]
      )
    )
  }

  # calculate cause specific hazards
  for (i in seq_rhs) {
    ped[[paste0("hazard", i)]] <-  exp(eval(F_rhs[[i]], ped))
  }
  ped[["rate"]] <- reduce(ped[paste0("hazard", seq_rhs)], `+`)

  # simulate survival times
  sim_df <- ped %>%
    group_by(id) %>%
    mutate(
      time_pexp   = rpexp_msm(rate = .data$rate, t = .data$tstart),
      time_event  = ifelse(is.null(round), .data$time_pexp, round(.data$time_pexp, digits = round)),
      status      = 1L * (.data$time_event <= max(cut)),
      time        = pmin(.data$time_event, max(cut))
    ) %>%
    ungroup() %>%
    filter(.data$tstart < .data$time & .data$time <= .data$tend) # only a single row per id

  # Ziehe aus den möglichen hazards eins mit den entsprechenden Wahrscheinlichkeiten
  sim_df$type <- apply(
    sim_df[paste0("hazard", seq_rhs)], 1,
    function(probs) sample(seq_rhs, 1, prob = probs)
  )

  # Map single 'from' argument and sampled 'to' to columns
  sim_df <- sim_df %>%
    mutate(
      from = from,
      to   = to[type]
    )

  # Update clocks: t_<from>_<to> for each possible to
  for(to_j in to) {
    col <- paste0("t_", from, "_", to_j)
    sim_df[[col]] <- ifelse(
      sim_df$to == to_j,
      sim_df$time - sim_df$t,
      NA_real_
    )
  }

  sim_df %>%
    select(-one_of(c("tstart", "tend", "interval", "offset", "ped_status", "rate")))

}


#' Simulate data for multi-state modeling scenario
#' via consecutive calls to \code{sim_pexp_cr_msm}
#' where start time and state are updated for each individual
#'
#' @keywords internal
sim_pexp_msm <- function(formulas_list, data, cut, terminal_states, round = NULL, add_counterfactuals = TRUE) {

  max_time <- max(cut)

  # Helper function to extract the "from", "to", and "formula" from each list element.
  extract_transition <- function(x) {
    from_val <- if (!is.null(x$from)) x$from else x[[1]]
    to_val   <- if (!is.null(x$to)) x$to else x[[2]]
    form_val <- if (!is.null(x$formula)) x$formula else x[[3]]
    list(from = from_val, to = to_val, formula = form_val)
  }

  # Extract all transitions
  trans_list  <- lapply(formulas_list, extract_transition)
  transitions <- vapply(trans_list, function(tr) paste0(tr$from, "_", tr$to), character(1))
  from_states <- unique(vapply(trans_list, function(tr) as.character(tr$from), character(1)))

  # Add one column per transition: time_<from><to>
  for(tr in transitions) {
    data[[paste0("t_", tr)]] <- NA_real_
  }
  # Add per-state clocks (excluding initial state 0):
  for(fs in from_states) {
    if (fs != "0") {
      data[[paste0("t_until_", fs)]] <- NA_real_
      data[[paste0("t_",          fs)]] <- NA_real_
    }
  }

  # Initialize an empty list to accumulate simulated transitions.
  sim_results <- list()

  # current_data holds individuals still at risk for further transitions.
  current_data <- data

  while (nrow(current_data) > 0) {

    # For each unique state in the current data:
    current_states <- as.character(unique(current_data$from))

    # For accumulating data for the next round:
    next_round <- list()

    for (state in current_states) {
      # Select individuals currently in this state.
      state_data <- current_data[current_data$from == state, ]

      # Identify transitions possible from this state.
      possible_transitions <- Filter(function(x) {
        trans <- extract_transition(x)
        as.character(trans$from) == as.character(state)
      }, formulas_list)

      to <- unique(sapply(possible_transitions, function(x) {
        trans <- extract_transition(x)
        trans$to
      }))

      # Skip if no transitions are defined for this state.
      if (length(possible_transitions) == 0) next

      # Combine the hazard formulas for all transitions from this state into one formula.
      combined_rhs <- vapply(possible_transitions, function(x) {
        trans <- extract_transition(x)
        deparse1(trans$formula[[2]])
      }, character(1))
      combined_formula <- as.formula(paste("~", paste(combined_rhs, collapse = " | ")))

      # Simulate waiting times for all competing transitions from this state.
      sim_df <- sim_pexp_cr_msm(combined_formula, state_data, state, to, cut, round = round)

      # pick up the t_until_* columns that were on state_data
      t_until_cols <- grep("^t_until_", names(state_data), value = TRUE)

      if (length(t_until_cols) > 0) {
        sim_df <- sim_df %>%
          # bring in the old columns, suffixed ".old"
          left_join(
            state_data %>% select(id, all_of(t_until_cols)),
            by     = "id",
            suffix = c("", ".old")
          )

        # for each t_until_<k>, if the freshly simulated value is NA, replace by the old one
        for (col in t_until_cols) {
          oldcol <- paste0(col, ".old")
          sim_df[[col]] <- ifelse(
            is.na(sim_df[[col]]),
            sim_df[[oldcol]],
            sim_df[[col]]
          )
        }

        # drop the helper ".old" columns
        sim_df <- sim_df %>% select(-ends_with(".old"))
      }

      # Combine the result from this state with the overall results.
      sim_results[[as.character(state)]] <- sim_df

      # Identify individuals who experienced an event (status==1),
      # who are not censored (time < max_time), and whose new state is not terminal.
      successful <- sim_df %>%
        filter(status == 1, time < max_time, !(to %in% terminal_states))

      # Update these individuals for the next round:
      # The new starting time is the simulated event time (in column 'time'),
      # and their current state is updated to the 'to' state.
      if (nrow(successful) > 0) {
        # pull everything you need from sim_df itself
        updated <- sim_df %>%
          filter(id %in% successful$id) %>%
          transmute(
            id,
            x1,                             # or any covariates you need
            from = to,                      # new state
            t    = time,                    # clock reset
            # carry forward all t_until_* columns
            across(starts_with("t_until_")),
            # carry forward all per‐transition clocks
            across(matches("^t_[0-9]+_[0-9]+"))
          )
        next_round[[as.character(state)]] <- updated
      }
    }

    # Prepare for the next iteration:
    current_data <- bind_rows(next_round)
    if (nrow(current_data) == 0) break
  }

  # Combine all simulated transitions into the final multi-state data frame.
  final_sim_df <- bind_rows(sim_results) %>%
    select(-one_of(c("type", "time_pexp", "time_event", "t_0", "t_until_0")),
           -dplyr::contains("hazard"),
           -starts_with("t_")) %>%
    mutate(transition = factor(paste0(from, "->", to))) %>%
    relocate(to, .after = from) %>%
    rename(tstart = t, tstop = time) %>%
    relocate(tstop, .after = tstart) %>%
    arrange(id, tstart)

  # Add counterfactual transitions if requested.
  if (add_counterfactuals) {
    final_sim_df <- final_sim_df %>% add_counterfactual_transitions()
  }

  return(final_sim_df)
}


#' Add censoring on top of simulated data
#' TBD: let user pass vector of times instead of imposing (parametric) censoring distribution
#' TBD: implement pexp for right censoring with covariate-dependent hazard
#' TBD: implement left-censoring
#' @keywords internal
add_censoring <- function(
  data,
  type = c("right", "interval", "left"),
  distribution = c("weibull", "exponential", "lognormal", "uniform"),
  parameters = NULL,
  round = NULL) {

  type <- match.arg(type)
  distribution <- match.arg(distribution)

  if (type == "right") {
    # Check parameters length based on distribution
    if (distribution %in% c("weibull", "lognormal")) {
      if (length(parameters) != 2) {
        stop("For 'weibull' or 'lognormal' distribution, 'parameters' must be of length 2 (e.g., shape & scale, or meanlog & sdlog).")
      }
    } else if (distribution == "exponential") {
      if (length(parameters) != 1) {
        stop("For 'exponential' distribution, 'parameters' must be of length 1 (i.e., rate).")
      }
    } else {
      stop("Unsupported distribution. Choose 'weibull', 'exponential', or 'lognormal'.")
    }

    data_right <- data %>%
      group_by(id) %>%
      # Draw one censoring time per individual
      mutate(censoring_time = case_when(
        distribution == "weibull"  ~ rweibull(1, shape = parameters[1], scale = parameters[2]),
        distribution == "exponential" ~ rexp(1, rate = parameters[1]),
        distribution == "lognormal" ~ rlnorm(1, meanlog = parameters[1], sdlog = parameters[2])
      )) %>%
      mutate(censoring_time = ifelse(is.null(round), censoring_time, round(censoring_time, digits = round))) %>%
      # Now, within the same grouping, order and determine censoring
      arrange(tstart) %>%
      mutate(row_index = row_number(),
            censor_flag = (tstart <= censoring_time & censoring_time < tstop),
            first_censor = ifelse(any(censor_flag), min(row_index[censor_flag]), Inf)) %>%
      filter(row_index <= first_censor) %>%
      mutate(tstop = ifelse(censor_flag, censoring_time, tstop),
            status = ifelse(censor_flag, 0, status)) %>%
      ungroup() %>%
      filter(tstart < tstop) %>%  # remove any row where tstart == tstop
      select(-c(censoring_time, censor_flag, row_index, first_censor)) %>%
      arrange(id, tstart)

    return(data_right)

  } else if (type == "left") {
    # Left-censoring branch applies to individuals with >1 transition.
    if (distribution %in% c("weibull", "lognormal")) {
      if (length(parameters) != 2) {
        stop("For 'weibull' or 'lognormal' distribution, 'parameters' must be of length 2 (e.g., shape & scale, or meanlog & sdlog).")
      }
    } else if (distribution == "exponential") {
      if (length(parameters) != 1) {
        stop("For 'exponential' distribution, 'parameters' must be of length 1 (i.e., rate).")
      }
    } else {
      stop("Unsupported distribution. Choose 'weibull', 'exponential', or 'lognormal'.")
    }

    # Process every individual via group_modify.
    # For groups with only one row, return unchanged.
    # For groups with >1 row, draw a left-censoring time L and adjust:
    # - If L <= first tstart: leave unchanged.
    # - Otherwise, find the first row i with L <= tstop[i],
    #   drop rows 1:(i-1) and update row i's tstart to L.
    # - If L > all tstop values, drop the individual.
    data_left <- data %>%
      group_by(id) %>%
      arrange(tstart) %>%
      group_modify(~ {
        df <- .x[order(.x$tstart), ]
        if(nrow(df) < 2) return(df)  # Only one row, leave unchanged.
        # Draw left-censoring time L for this individual.
        L <- if (distribution == "weibull") {
          rweibull(1, shape = parameters[1], scale = parameters[2])
        } else if (distribution == "exponential") {
          rexp(1, rate = parameters[1])
        } else if (distribution == "lognormal") {
          rlnorm(1, meanlog = parameters[1], sdlog = parameters[2])
        }
        if (!is.null(round)) {
          L <- round(L, digits = round)
        }
        # If L is less than or equal to the first observed tstart, leave unchanged.
        if (L <= df$tstart[1]) {
          return(df)
        } else {
          # Find the first row i such that L <= tstop[i]
          i <- which(L <= df$tstop)[1]
          if (is.na(i)) {
            # L exceeds all tstop values: drop the individual.
            return(tibble())
          } else {
            # Drop rows 1 to (i-1) and update row i's tstart to L.
            df_new <- df[i:nrow(df), ]
            df_new$tstart[1] <- L
            return(df_new)
          }
        }
      }) %>% ungroup()

    return(data_left)

  } else if (type == "interval") {
    if (distribution != "uniform") {
      stop("For interval censoring, only the 'uniform' distribution is supported.")
    }

    # Interval-censoring branch: for each individual with multiple transitions,
    # leave the first tstart and the final tstop unchanged.
    # For each intermediate interval j (1 <= j < k), draw an imputed event time from
    # [tstop[j], tstop[j+1]) and update:
    #   - row j's tstop becomes the imputed value,
    #   - row (j+1)'s tstart becomes the imputed value.
    data_interval <- data %>%
      arrange(id, tstart) %>%
      group_by(id) %>%
      group_modify(~ {
        df <- .x[order(.x$tstart), ]
        k <- nrow(df)
        if (k < 2) return(df)  # Nothing to impute if only one row.

        imputed <- numeric(k - 1)
        for (j in 1:(k - 1)) {
          if (is.null(round)) {
            imputed[j] <- runif(1, min = df$tstop[j], max = df$tstop[j + 1])
          } else {
            step <- 10^(-round)
            lower <- ceiling(df$tstop[j] * 10^round) / 10^round
            grid_max <- df$tstop[j + 1] - step
            if (grid_max < lower) {
              imputed[j] <- round(runif(1, min = df$tstop[j], max = df$tstop[j + 1]), round)
            } else {
              possible_vals <- seq(lower, grid_max, by = step)
              if (length(possible_vals) == 0) {
                imputed[j] <- round(runif(1, min = df$tstop[j], max = df$tstop[j + 1]), round)
              } else {
                imputed[j] <- sample(possible_vals, 1)
              }
            }
          }
        }
        df$tstop[1] <- imputed[1]
        if (k > 2) {
          for (j in 2:(k - 1)) {
            df$tstart[j] <- imputed[j - 1]
            df$tstop[j] <- imputed[j]
          }
        }
        df$tstart[k] <- imputed[k - 1]
        return(df)
      }) %>% ungroup()

    return(data_interval)

  } else {
    stop("Censoring type must be one of 'right', 'left', or 'interval'.")
  }
}


add_timeScales <- function(ped) {

    saved_attrs <- attributes(ped)

    out <- ped %>%
        # For each individual, determine the onset and progression ages:
        group_by(id) %>%
        mutate(
            t_until_1 =
                if(any(from == 1)) {
                    first(tstart[from == 1])
                } else if(any(to == 1)){
                    last(tend[to == 1])
                } else if(any(from== 2)){ # for those whose first transition is 2->
                    first(tstart[from == 2])
                } else {
                    0
                },
            t_until_2 =
                if(any(from == 2)) {
                    first(tstart[from == 2])
                } else if(any(to == 2)){
                    last(tend[to == 2])
                } else {
                    0
                },
        ) %>%
        ungroup() %>%
        # Now, compute the onset columns,
        # forcing them to 0 for rows with 0-> transitions.
        mutate(
            t_until_1 = case_when(
                from == 0 ~ 0,
                TRUE ~ t_until_1
            ),
            tstart_1 = case_when(
                t_until_1 == 0 ~ 0,
                from == 0 ~ 0,
                TRUE ~ tstart - t_until_1
            ),
            t_1 = case_when(
                t_until_1 == 0 ~ 0,
                from == 0 ~ 0,
                TRUE ~ tend - t_until_1
            ),
            # And compute the progression columns:
            # For individuals with no 2-> transition (i.e. t_until_2 == 0)
            # or for rows that occur before the progression event,
            # we set progression values to 0.
            t_until_2 = case_when(
                from %in% c(0, 1) ~ 0,
                TRUE ~ t_until_2
            ),
            tstart_2 = case_when(
                t_until_2 == 0 ~ 0,
                from %in% c(0, 1) ~ 0,
                TRUE ~ tstart - t_until_2
            ),
            t_2 = case_when(
                t_until_2 == 0 ~ 0,
                from %in% c(0, 1) ~ 0,
                TRUE ~ tend - t_until_2
            ),
            t_1_2 = case_when(
                from == 2 & t_until_2 > 0 & t_until_1 > 0 ~ t_until_2 - t_until_1,
                TRUE ~ 0
            )
        )

    # Reapply all original attributes (excluding names, row.names, and class which mutate preserves)
    attr_names <- setdiff(names(saved_attrs), c("names", "row.names", "class"))
    for (a in attr_names) {
    attr(out, a) <- saved_attrs[[a]]
    }

    return(out)
}


# add_transVars <- function(ped) {
#   out <- ped %>%
#     mutate(
#       trans_to_4 = factor(ifelse(to == 4, 1, 0), levels = c(0, 1)),
#       trans_after_1 = factor(
#         case_when(
#           from == 0 ~ "none",
#           transition %in% c("1->2", "2->3") ~ "progression",
#           transition %in% c("1->4", "2->4") ~ "4"
#         ),
#           levels = c("none", "progression", "4")
#       ),
#       trans_after_2 = factor(
#         case_when(
#           from %in% c(0, 1) ~ "none",
#           transition == "2->3" ~ "3",
#           transition == "2->4" ~ "4"
#         ),
#         levels = c("none", "3", "4")
#       ),
#       trans_after_1_strat = factor(
#         case_when(
#           from == 0 ~ "none",
#           TRUE ~ transition
#         ),
#         levels = c("none", "1->2", "2->3", "1->4", "2->4")
#       )
#     )

#   return(out)
# }


add_transVars <- function(ped) {
  out <- ped %>%
    mutate(
      trans_to_3 = factor(ifelse(to == 3, 1, 0), levels = c(0, 1)),
      trans_after_1 = factor(
        case_when(
          from == 0 ~ "none",
          TRUE ~ transition
        ),
        levels = c("none", "1->2", "1->3")
      )
    )

  return(out)
}
