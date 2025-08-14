library(data.table)
library(dplyr)
library(dtplyr)
library(tidyr)
library(rlang) # for sim_pexp
library(Formula) # for sim_pexp
library(lazyeval) # for sim_pexp
library(purrr) # for sim_pexp
library(stringr)
library(tibble)
library(survival)
library(mgcv)
library(pammtools)
library(flexsurv)
library(ggplot2)
library(scales)

# data prep ----
wrapper_sim <- function(
  data,
  job,
  formulas_dgp,
  terminal_states,
  cut,
  n = 2500,
  round = 2,
  cens_type = c("none", "right", "interval", "left"),
  cens_dist = c("weibull", "exponential", "lognormal", "uniform"),
  cens_params = NULL){

  df <- data.frame(
    id          = seq_len(n),
    from        = 0L,
    t           = 0
  )

  vars <- get_x_vars(formulas_dgp)
  for (var in vars) {
    df[[var]] <- rbinom(n, 1, 0.5)
  }

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
    ped = ped
  )

  return(out)
}


wrapper_fe <- function(
  data,
  job,
  instance,
  formula,
  bs_0 = "ps",
  bs = "ps",
  k = 10) {

  formulas_dgp <- data$formulas_dgp
  ped <- instance$ped

  if(bs != "fs") {
    ped$trans_after_1 <- factor(ped$trans_after_1,
                                levels = c("none", "1->2", "1->3"),
                                ordered = TRUE)
  }

  formula <- convert_formula(formula, bs_0, bs, k)

  # Extract transitions from the formulas_dgp
  transitions <- map_chr(formulas_dgp, ~ {
    from <- .x$from
    to   <- .x$to
    paste0(from, "->", to)
  })

  # Extract all beta_1 coefficients from formulas_dgp
  beta_1_vars <- formulas_dgp %>%
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

  # Fit the model
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
      summary_mod$p.coeff["x1"] + summary_mod$p.coeff[paste0("transition", transition_name, ":", "x1")]
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

wrapper_bh <- function(
  data,
  job,
  instance,
  formula,
  bs_0 = "ps",
  bs = "ps",
  k = 10,
  ci = TRUE) {

  formulas_dgp <- data$formulas_dgp
  ped          <- instance$ped
  cut          <- data$cut

  if(bs != "fs") {
    ped$trans_after_1 <- factor(ped$trans_after_1,
                                levels = c("none", "1->2", "1->3"),
                                ordered = TRUE)
  }

  formula <- convert_formula(formula, bs_0, bs, k)

  # Extract transitions from the formulas_dgp
  transitions <- map_chr(formulas_dgp, ~ {
    from <- .x$from
    to   <- .x$to
    paste0(from, "->", to)
  })

  # Fit the model
  mod <- bam(as.formula(formula)
                      , data = ped
                      , family = poisson()
                      , offset = offset
                      , discrete = T
                      , method = "fREML"
  )

  dgp_tbl <- tibble(
    from = as.character(map_int(formulas_dgp, "from")),
    to   = as.character(map_int(formulas_dgp, "to")),
    expr = map(formulas_dgp, "formula")        # each is a one‐sided formula: ~ f_0(t)+…
  )

  nd <- make_newped(ped, cut, mod, bs, ci) %>%
    left_join(dgp_tbl, by = c("from", "to")) %>%
    rowwise() %>%
    mutate(
      loghazard_true = {
        rhs <- expr[[2]] # strip off the “~” and evaluate the RHS in the current row
        eval(rhs, envir = list(
          t          = tend,
          t_1        = t_1,
          t_until_1  = t_until_1
        ))
      },
      hazard_true = exp(loghazard_true)
    ) %>%
    ungroup()

  # split by 'from'
  nd0 <- nd %>% filter(from == "0")
  nd1 <- nd %>% filter(from != "0")

  # process from=0: group only by transition, arrange without t_until_1
  nd0 <- nd0 %>%
    arrange(transition, tend) %>%
    group_by(transition) %>%
    rename(
      cumu_hazard_emp = cumu_hazard,
      trans_prob_emp  = trans_prob
    ) %>%
    mutate(
      cumu_hazard = cumsum(hazard_true * intlen)
    ) %>%
    add_trans_prob(mod, ci = FALSE) %>%
    rename(
      cumu_hazard_true = cumu_hazard,
      trans_prob_true  = trans_prob,
      cumu_hazard      = cumu_hazard_emp,
      trans_prob       = trans_prob_emp
    ) %>%
    ungroup()

  # process from=1: keep the original grouping by (t_until_1, transition)
  nd1 <- nd1 %>%
    arrange(transition, t_until_1, tend) %>%
    group_by(transition, t_until_1) %>%
    rename(
      cumu_hazard_emp = cumu_hazard,
      trans_prob_emp  = trans_prob
    ) %>%
    mutate(
      cumu_hazard = cumsum(hazard_true * intlen)
    ) %>%
    add_trans_prob(mod, ci = FALSE) %>%
    rename(
      cumu_hazard_true = cumu_hazard,
      trans_prob_true  = trans_prob,
      cumu_hazard      = cumu_hazard_emp,
      trans_prob       = trans_prob_emp
    ) %>%
    ungroup()

  # bind them back
  nd <- bind_rows(nd0, nd1)

  if(ci) {
    nd <- nd %>%
      rename(
        cumu_hazard_lower = cumu_lower,
        cumu_hazard_upper = cumu_upper,
        trans_prob_lower = trans_lower,
        trans_prob_upper = trans_upper
      ) %>%
      mutate(
        loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)),
        hazard_cov = as.integer((hazard_true >= hazard_lower) & (hazard_true <= hazard_upper)),
        cumu_hazard_cov = as.integer((cumu_hazard_true >= cumu_hazard_lower) & (cumu_hazard_true <= cumu_hazard_upper)),
        trans_prob_cov =  as.integer((trans_prob_true >= trans_prob_lower) & (trans_prob_true <= trans_prob_upper))
      )
  }

  if(bs != "fs") {
    nd$trans_after_1 <- factor(nd$trans_after_1,
                                levels = c("none", "1->2", "1->3"),
                                ordered = FALSE)
  }

  return(nd)

}


wrapper_bh_terms <- function(
  data,
  job,
  instance,
  formula) {

  formulas_dgp <- data$formulas_dgp
  ped          <- instance$ped
  cut          <- data$cut

  if(bs != "fs") {
    ped$trans_after_1 <- factor(ped$trans_after_1,
                                levels = c("none", "1->2", "1->3"),
                                ordered = TRUE)
  }

  formula <- convert_formula(formula, bs_0, bs, k)

  # Extract transitions from the formulas_dgp
  transitions <- map_chr(formulas_dgp, ~ {
    from <- .x$from
    to   <- .x$to
    paste0(from, "->", to)
  })

  # Fit the model
  mod <- bam(as.formula(formula)
                      , data = ped
                      , family=poisson()
                      , offset=offset
                      , discrete = T
                      , method = "fREML"
  )

png(file.path("loghazards_from_mod.png"), width = 1000, height = 800)
plot(mod, pages = 1)
dev.off()

  dgp_tbl <- tibble(
    from = map_chr(formulas_dgp, ~ as.character(.x$from)),
    to   = map_chr(formulas_dgp, ~ as.character(.x$to)),
    expr = map(formulas_dgp, "formula")   # each is ~ f_*(t…) + …
  ) %>%
    mutate(
      terms = map(expr, ~ {
        # pull out all the names in the RHS
        vars <- all.vars(.x)
        # drop any that are:
        #  - the “until” time
        #  - your f_*/g_* functions
        #  - the beta constants
        keep <- vars[! str_detect(vars, "_until_")
                    & ! str_detect(vars, "^(f_|g_)")
                    & ! str_detect(vars, "^beta")]
        unique(keep)
      })
    )

  nd <- make_newped_terms(ped, cut, mod, dgp_tbl, from_state = 1, ci) %>%
    mutate(
      loghazard_true = case_when(
        transition == "1->2" ~ f_1(t_1) + beta_0_12,
        transition == "1->3" ~ g_1(t_1) + beta_0_13,
        # transition == "1->2" ~ f_1(t_1),
        # transition == "1->3" ~ g_1(t_1),
        TRUE                 ~ NA_real_
      )
  )
  # # group by transition and subtract off the per‐transition mean
  # group_by(transition) %>%
  # mutate(
  #   true_mean        = mean(loghazard_true_uncentered, na.rm = TRUE),
  #   loghazard_true = loghazard_true_uncentered - true_mean
  # ) %>%
  # ungroup()

  if(ci) {
    nd <- nd %>%
      mutate(loghazard_cov = as.integer((loghazard_true >= loghazard_lower) & (loghazard_true <= loghazard_upper)))
  }

trans <- "1->2"

df <- nd %>% filter(transition == trans) %>% group_by(t_1) %>% slice(1)

  gg <- ggplot(df, aes(x = t_1, y = loghazard)) +
    geom_line(color = "black") +
    geom_ribbon(aes(ymin = loghazard_lower, ymax = loghazard_upper), fill = "black", alpha = 0.2) +
    geom_line(aes(y = loghazard_true), color = "orange", linetype = "dashed") +
    xlab("time") +
    ylab("loghazard")

  ggsave(
    filename = file.path("loghazard.png"),
    plot     = gg,
    width    = 10, height = 8
  )

  return(nd)

}


convert_formula <- function(formula, bs_0, bs, k) {
  # Step i: Convert formula to character
  f_txt <- paste(deparse(formula), collapse = " ")

  # Step ii: Remove "by = " before trans_after_1 if bs == "fs"
  if (bs == "fs") {
    f_txt <- gsub("\\bby\\s*=\\s*trans_after_1\\b", "trans_after_1", f_txt)
  }

  # Step iii: Evaluate and bake in the values of bs_0, bs, and k
  # Replace bs_0 with its evaluated value (only when it's a standalone variable, not parameter name)
  f_txt <- gsub("\\bbs\\s*=\\s*bs_0\\b", paste0('bs = "', bs_0, '"'), f_txt)

  # Replace bs with its evaluated value (only when it's a value, not parameter name)
  f_txt <- gsub("\\bbs\\s*=\\s*bs\\b", paste0('bs = "', bs, '"'), f_txt)

  # Replace k with its evaluated value (only when it's a value, not parameter name)
  f_txt <- gsub("\\bk\\s*=\\s*k\\b", paste0('k = ', k), f_txt)

  # Step iv: Convert back to formula
  new_formula <- as.formula(f_txt)

  return(new_formula)
}


make_newped <- function(ped, cut, mod, bs, ci = FALSE) {

  from_states <- sort(unique(ped$from))
  transitions <- sort(unique(ped$transition))

  ped_new_in <- expand_grid(
    transition = transitions,
    tend = sort(cut)[-1], # all but the first cut
    t_until_1 = sort(cut)[-length(cut)] # all but the last cut
  )

  ped_new_in <- ped_new_in %>%
    mutate(
      from = as.integer(sub("([0-9])->.*","\\1", transition)),
      to   = as.integer(sub(".*->([0-9])","\\1", transition)),
      tstart = tend - diff(cut)[1], # using the difference between consecutive values of cut
      intlen = tend - tstart,
      t_1 = ifelse(from == 0, 0, tend - t_until_1),
      t_until_1 = ifelse(from == 0, 0, t_until_1)
    ) %>%
    add_transVars()

  if(bs != "fs") {
    ped_new_in$trans_after_1 <- factor(ped_new_in$trans_after_1,
                                  levels = c("none", "1->2", "1->3"),
                                  ordered = TRUE)
  }

  ped_new <- bind_rows(
    lapply(from_states, function(from) {
      ped_new_trans_in <- ped_new_in %>% filter(from == !!from)

      if (from == "0") {
        ped_new_trans <- ped_new_trans_in %>%
          distinct() # for transition 0->1, for a given tend, all rows are identical (because t_until_1 = 0 always); same for 0->4
      } else if (from == "1") {
        ped_new_trans <- ped_new_trans_in %>%
          distinct() %>%
          filter(tstart >= t_until_1) %>%
          filter(t_until_1 < 9.9)
      }

      return(ped_new_trans)
    })
  )

  # split off the two sets of transitions
  ped_from0 <- ped_new %>% filter(from == 0)
  ped_from1 <- ped_new %>% filter(from != 0)

  # process the from‐0 transitions: group/arrange only by transition & tend
  ped_from0 <- ped_from0 %>%
    group_by(transition) %>%
    arrange(transition, tend) %>%
    add_hazard(mod, type = "link", ci = ci) %>%
    rename(loghazard = hazard) %>%
    { if(ci) rename(., loghazard_lower = ci_lower, loghazard_upper = ci_upper) else . } %>%
    add_hazard(mod, ci = ci) %>%
    { if(ci) rename(., hazard_lower = ci_lower, hazard_upper = ci_upper) else . } %>%
    add_cumu_hazard(mod, ci = ci) %>%
    add_trans_prob(mod, ci = ci)

  # process the from‐1 transitions: keep t_until_1 in the grouping
  ped_from1 <- ped_from1 %>%
    group_by(t_until_1, transition) %>%
    arrange(t_until_1, transition, tend) %>%
    add_hazard(mod, type = "link", ci = ci) %>%
    rename(loghazard = hazard) %>%
    { if(ci) rename(., loghazard_lower = ci_lower, loghazard_upper = ci_upper) else . } %>%
    add_hazard(mod, ci = ci) %>%
    { if(ci) rename(., hazard_lower = ci_lower, hazard_upper = ci_upper) else . } %>%
    add_cumu_hazard(mod, ci = ci) %>%
    add_trans_prob(mod, ci = ci)

  # put them back together, then ungroup + factor / from / to recoding
  ped_new <- bind_rows(ped_from0, ped_from1) %>%
    ungroup() %>%
    mutate(
      transition = factor(transition, levels = c("0->1","1->2","0->3","1->3")),
      from       = sub("([0-9])->.*","\\1", transition),
      to         = sub(".*->([0-9])","\\1", transition)
    )

  # ped_new <- ped_new %>%
  #   group_by(t_until_1, transition) %>%
  #   arrange(t_until_1, transition, tend) %>%
  #   add_hazard(mod, type = "link", ci = ci) %>%
  #   rename(loghazard = hazard)

  # if(ci) {
  #   ped_new <- ped_new %>%
  #     rename(
  #       loghazard_lower = ci_lower,
  #       loghazard_upper = ci_upper
  #     )
  # }

  # ped_new <- ped_new %>%
  #   add_hazard(mod, ci = ci)

  # if(ci) {
  #   ped_new <- ped_new %>%
  #     rename(
  #       hazard_lower = ci_lower,
  #       hazard_upper = ci_upper
  #     )
  # }

  # ped_new <- ped_new %>%
  #   add_cumu_hazard(mod, ci = ci) %>%
  #   add_trans_prob(mod, ci = ci) %>%
  #   ungroup() %>%
  #   mutate(
  #     transition = factor(transition, levels = c("0->1", "1->2", "0->3", "1->3")),
  #     from = as.character(sub("([0-9])->.*","\\1", transition)),
  #     to   = as.character(sub(".*->([0-9])","\\1", transition))
  #   )

  return(ped_new)
}


# make_newped_terms <- function(ped, cut, mod, dgp_tbl, from_state, ci = FALSE) {

#   from_states <- sort(unique(ped$from))
#   transitions <- sort(unique(ped$transition))

#   ped_new_in <- expand_grid(
#     transition = transitions,
#     tend = sort(cut)[-1], # all but the first cut
#     t_until_1 = sort(cut)[-length(cut)] # all but the last cut
#   )

#   ped_new_in <- ped_new_in %>%
#     mutate(
#       from = as.integer(sub("([0-9])->.*","\\1", transition)),
#       to   = as.integer(sub(".*->([0-9])","\\1", transition)),
#       tstart = tend - diff(cut)[1], # using the difference between consecutive values of cut
#       intlen = tend - tstart,
#       t_until_1 = ifelse(from == 0, 0, t_until_1),
#       t_1 = ifelse(from == 0, 0, tend - t_until_1)
#     ) %>%
#     add_transVars()

#   ped_new <- bind_rows(
#     lapply(from_states, function(f) {
#       df <- ped_new_in %>% filter(from == f) %>% distinct()
#       if (f == 1) {
#         df <- df %>%
#           filter(tstart >= t_until_1, t_until_1 < max(cut))
#       }
#       df
#     })
#   )

#   ped_new <- ped_new %>% filter(from == from_state)

#   ped_new <- ped_new %>%
#     mutate(
#       from = as.character(from),
#       to   = as.character(to)
#     ) %>%
#     left_join(dgp_tbl, by = c("from","to"))

#   valid_tr <- dgp_tbl %>%
#     # only those transitions that have a t_1 term
#     filter(map_lgl(terms, ~ "t_1" %in% .x)) %>%
#     transmute(
#       tr          = paste0(from, "->", to),
#       par_term    = paste0("transition", tr),
#       smooth_term = "t_1"
#     )

#   ped_new_terms <- purrr::map_dfr(valid_tr$tr, function(tr) {
#     # subset to this transition
#     df <- ped_new %>% filter(as.character(transition) == tr)

#     # build the vector of three term‐names
#     tt <- valid_tr %>% filter(tr == !!tr)
#     term_vec <- c("(Intercept)", tt$par_term, tt$smooth_term)
#     # term_vec <- tt$smooth_term

#     # one call to add_term() for the sum of all three
#     df %>%
#       add_term(
#         object        = mod,
#         term          = term_vec,
#         se            = ci,
#         unconditional = TRUE
#       ) %>%
#       # rename the returned 'fit' into your log‐hazard
#       rename(loghazard = fit) %>%
#       # rename CIs only if requested
#       { if (ci) rename(.,
#           loghazard_lower = ci_lower,
#           loghazard_upper = ci_upper
#         ) else . }
#   })


#   ped_new_terms <- ped_new_terms %>%
#     select(-c(expr, terms))

#   return(ped_new_terms)
# }


create_fe_table <- function(df, grouping_vars = character()) {
  # 0) Validate grouping_vars is a character vector
  if (!is.character(grouping_vars)) {
    stop("`grouping_vars` must be a character vector")
  }

  # 1) Check that each grouping_var actually exists in df
  missing_vars <- setdiff(grouping_vars, colnames(df))
  if (length(missing_vars) > 0) {
    stop("The following `grouping_vars` are not columns in `df`: ",
         paste(missing_vars, collapse = ", "))
  }

  # 2) Build the full grouping key: always include `variable`
  groups <- c("variable", grouping_vars)

  # 3) For each group, check true_value is constant
  inconsistent <- df %>%
    group_by(across(all_of(groups))) %>%
    summarize(n_true = n_distinct(true_value), .groups = "drop") %>%
    filter(n_true > 1)

  if (nrow(inconsistent) > 0) {
    bad <- inconsistent %>%
      unite("grp", all_of(groups), sep = " | ") %>%
      pull(grp)
    stop("Inconsistent `true_value` in groups:\n  ",
         paste(bad, collapse = "\n  "))
  }

  # 4) Summarize
  df %>%
    group_by(across(all_of(groups))) %>%
    summarize(
      true_value       = first(true_value),
      average_estimate = mean(estimate),
      average_coverage = mean(coverage),
      bias             = average_estimate - true_value,
      .groups = "drop"
    ) %>%
    # ensure the columns come in the desired order:
    select(variable, all_of(grouping_vars),
           true_value, average_estimate, average_coverage, bias)
}


create_coverage_summary_fe <- function(df, grouping_vars = character()) {
  # 0) Validate grouping_vars is a character vector
  if (!is.character(grouping_vars)) {
    stop("`grouping_vars` must be a character vector")
  }

  # 1) Check that each grouping_var actually exists in df
  missing_vars <- setdiff(grouping_vars, colnames(df))
  if (length(missing_vars) > 0) {
    stop("The following `grouping_vars` are not columns in `df`: ",
         paste(missing_vars, collapse = ", "))
  }

  # 2) Build the full grouping key: always include `problem` and `variable`
  # --- THIS LINE IS CORRECTED ---
  groups <- c("problem", "variable", grouping_vars)

  # 3) For each group, check true_value is constant
  inconsistent <- df %>%
    group_by(across(all_of(groups))) %>%
    summarize(n_true = n_distinct(true_value), .groups = "drop") %>%
    filter(n_true > 1)

  if (nrow(inconsistent) > 0) {
    bad <- inconsistent %>%
      unite("grp", all_of(groups), sep = " | ") %>%
      pull(grp)
    stop("Inconsistent `true_value` in groups:\n  ",
         paste(bad, collapse = "\n  "))
  }

  # 4) Summarize
  summary_df <- df %>%
    group_by(across(all_of(groups))) %>%
    summarize(
      true_value       = first(true_value),
      average_estimate = mean(estimate),
      n_total          = dplyr::n(),
      n_covered        = sum(coverage, na.rm = TRUE),
      average_coverage = n_covered / n_total,
      bias             = average_estimate - true_value,
      .groups = "drop"
    )

  # 5) Calculate exact binomial confidence intervals
  final_summary_df <- summary_df %>%
    mutate(
      ci = purrr::map2(n_covered, n_total, ~ stats::binom.test(.x, .y)$conf.int),
      coverage_lower = purrr::map_dbl(ci, 1),
      coverage_upper = purrr::map_dbl(ci, 2)
    ) %>%
    # --- THIS LINE IS CORRECTED ---
    # Select final columns in the desired order, now including `problem`
    select(
      problem, variable, all_of(grouping_vars),
      true_value, average_estimate,
      coverage = average_coverage,
      coverage_lower,
      coverage_upper,
      bias
    )

  return(final_summary_df)
}


create_bh_table <- function(df,
                            grouping_vars = character(),
                            time_scales   = character()) {

  ## 1) ─ validation ------------------------------------------------------
  if (!is.character(grouping_vars))
    stop("`grouping_vars` must be a character vector.")
  miss_gv <- setdiff(grouping_vars, names(df))
  if (length(miss_gv))
    stop("Unknown grouping_vars: ", paste(miss_gv, collapse = ", "))

  if (!is.character(time_scales) || length(time_scales) < 1)
    stop("`time_scales` must be a non-empty character vector.")
  miss_ts <- setdiff(time_scales, names(df))
  if (length(miss_ts))
    stop("Unknown time_scales: ", paste(miss_ts, collapse = ", "))

  ## 2) ─ data.table setup ----------------------------------------------
  setDT(df)                                # in-place; no copy
  grp_keys <- c(grouping_vars, "transition", time_scales)

  ## 3) ─ summarise -------------------------------------------------------
  summaryDT <- df[ , {

      ## ------- constants & helpers
      N   <- .N
      x_lh <- sum(loghazard_cov)
      x_hz <- sum(hazard_cov)
      x_ch <- sum(cumu_hazard_cov)
      x_tp <- sum(trans_prob_cov)

      ci_lh <- stats::binom.test(x_lh, N, p = 0.5)$conf.int
      ci_hz <- stats::binom.test(x_hz, N, p = 0.5)$conf.int
      ci_ch <- stats::binom.test(x_ch, N, p = 0.5)$conf.int
      ci_tp <- stats::binom.test(x_tp, N, p = 0.5)$conf.int

      ## ------- means, bias, rmse as before -----------------------------
      t_lh <- first(loghazard_true);   m_lh <- mean(loghazard)
      t_hz <- first(hazard_true);      m_hz <- mean(hazard)
      t_ch <- first(cumu_hazard_true); m_ch <- mean(cumu_hazard)
      t_tp <- first(trans_prob_true);  m_tp <- mean(trans_prob)

      .(
        ## means
        loghazard_mean   = m_lh,
        hazard_mean      = m_hz,
        cumu_hazard_mean = m_ch,
        trans_prob_mean  = m_tp,

        ## true values
        loghazard_true   = t_lh,
        hazard_true      = t_hz,
        cumu_hazard_true = t_ch,
        trans_prob_true  = t_tp,

        ## bias
        loghazard_bias       = m_lh - t_lh,
        hazard_bias          = m_hz - t_hz,
        cumu_hazard_bias     = m_ch - t_ch,
        trans_prob_bias      = m_tp - t_tp,

        ## coverage proportions
        loghazard_coverage   = x_lh / N,
        hazard_coverage      = x_hz / N,
        cumu_hazard_coverage = x_ch / N,
        trans_prob_coverage  = x_tp / N,

        ## coverage CIs (exact binomial)
        loghazard_coverage_lower   = ci_lh[1],
        loghazard_coverage_upper   = ci_lh[2],
        hazard_coverage_lower      = ci_hz[1],
        hazard_coverage_upper      = ci_hz[2],
        cumu_hazard_coverage_lower = ci_ch[1],
        cumu_hazard_coverage_upper = ci_ch[2],
        trans_prob_coverage_lower  = ci_tp[1],
        trans_prob_coverage_upper  = ci_tp[2],

        ## RMSE
        loghazard_rmse       = sqrt(mean((loghazard    - t_lh)^2)),
        hazard_rmse          = sqrt(mean((hazard       - t_hz)^2)),
        cumu_hazard_rmse     = sqrt(mean((cumu_hazard  - t_ch)^2)),
        trans_prob_rmse      = sqrt(mean((trans_prob   - t_tp)^2))
      )
    }, by = grp_keys]

  summaryDT[]
}



create_bh_linePlots <- function(
  data,
  trans,
  grouping_vars = character(),
  scale = c("loghazard", "hazard", "cumu_hazard", "trans_prob"),
  font_size = 14,
  alpha = 1) {

  if (!is.character(grouping_vars))
    stop("`grouping_vars` must be a character vector.")

  match.arg(scale)

  y_col <- scale
  y_col_true <- paste0(scale, "_true")
  y_col_cov <- paste0(scale, "_cov")

  y_lab <- switch(
    scale,
    loghazard = "Log-hazard",
    hazard = "Hazard",
    cumu_hazard = "Cumulative Hazard",
    trans_prob = "Transition Probability"
  )

  required_cols <- c(
    grouping_vars,
    "job.id",
    "tend",
    scale,
    paste0(scale, "_true")
  )

  if (!all(required_cols %in% colnames(data))) {
    stop("Data is missing required columns.")
  }

  plot_data <- data %>%
    filter(transition == trans)

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

    rng       <- range(df_grp[[y_col]], df_grp[[y_col_true]], na.rm = TRUE)
    rng_min   <- rng[1]
    rng_span  <- diff(rng)
    ref_y     <- rng_min + 0.95 * rng_span

    cov_line  <- df_grp |>
      dplyr::group_by(tend) |>
      dplyr::summarise(
        cov_mean = mean(.data[[y_col_cov]], na.rm = TRUE),
        .groups  = "drop"
      ) |>
      dplyr::mutate(cov_y = rng_min + cov_mean * rng_span)

    p <- ggplot(df_grp, aes(x = tend, y = .data[[y_col]])) +
      geom_line(aes(group = job.id), alpha = alpha) +
      geom_line(aes(y = .data[[y_col_true]], colour = "truth"), linewidth = 1.5) +
      geom_smooth(aes(col = "average estimate"), method = "gam",
            formula = y ~ s(x), se = FALSE) +
      geom_line(data = cov_line,                 # dotted blue coverage curve
                aes(x = tend, y = cov_y),
                colour = "blue", linetype = "dotted") +
      geom_hline(yintercept = ref_y, linetype = "dashed", color = "red") +
      scale_y_continuous(
        name    = y_lab,
        sec.axis = sec_axis(~ (. - rng_min) / rng_span,
                name   = "Coverage",
                breaks = seq(0, 1, 0.2),
                labels = scales::percent_format(accuracy = 1))
      ) +
      scale_color_brewer("", palette = "Dark2") +
      xlab("Time") +
      ylab(y_lab) +
      # ggtitle(paste("Group:", grp)) +
      theme_bw() +
      theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      # legend.justification = c("left", "top"),
      legend.text = element_text(size = font_size),
      legend.title = element_text(size = font_size),
      axis.text = element_text(size = font_size),
      axis.title = element_text(size = font_size)
      )
  })

  names(plots) <- unique_groups
  return(plots)
}


create_bh_slicePlots <- function(
  data,
  trans,
  grouping_vars   = character(),
  scale           = c("loghazard", "hazard", "cumu_hazard", "trans_prob"),
  time_axis       = "tend",
  slice_axis      = "t_until_1",
  slice_points,
  ncol_facets     = 2,
  font_size       = 14,
  alpha           = 1
) {
  ## ---- checks ----------------------------------------------------------
  if (!is.character(grouping_vars))
    stop("`grouping_vars` must be a character vector.")
  scale <- match.arg(scale)
  if (missing(slice_points) || length(slice_points) == 0)
    stop("Please supply non-empty `slice_points`.")

  needed <- c(grouping_vars,
              "job.id", time_axis, slice_axis,
              scale, paste0(scale, "_true"), paste0(scale, "_cov"))
  miss <- setdiff(needed, names(data))
  if (length(miss))
    stop("`data` is missing: ", paste(miss, collapse = ", "))

  ## ---- data prep (all tidyverse) --------------------------------------
  plot_data <- data %>%
    filter(transition == trans) %>%
    { if (length(grouping_vars))
        unite(., grouping, all_of(grouping_vars), remove = FALSE)
      else
        mutate(., grouping = "All") } %>%
    filter(.data[[slice_axis]] %in% slice_points) %>%
    { if (slice_axis == time_axis)
        filter(., .data[[time_axis]] >= .data[[slice_axis]])
      else .
    } %>%
    mutate(
      slice_lab = factor(
        .data[[slice_axis]],
        levels = slice_points,
        labels = paste0(slice_axis, " = ", slice_points)
      )
    )

  ## column helpers
  y_col      <- scale
  y_col_true <- paste0(scale, "_true")
  y_col_cov  <- paste0(scale, "_cov")
  y_lab      <- scale

  ## ---- plot builder ----------------------------------------------------
  make_plot <- function(df) {
    rng_min  <- min(df[[y_col]], df[[y_col_true]], na.rm = TRUE)
    rng_max  <- max(df[[y_col]], df[[y_col_true]], na.rm = TRUE)
    rng_span <- rng_max - rng_min
    ref_y    <- rng_min + 0.95 * rng_span

    cov_line <- df %>%
      group_by(.data[[time_axis]], slice_lab) %>%
      summarise(cov_mean = mean(.data[[y_col_cov]], na.rm = TRUE),
                .groups  = "drop") %>%
      mutate(cov_y = rng_min + cov_mean * rng_span)

    ggplot(df, aes(x = .data[[time_axis]], y = .data[[y_col]])) +
      geom_line(aes(group = job.id),
                colour = "black", linewidth = 0.5, alpha = alpha) +
      geom_line(aes(y = .data[[y_col_true]], colour = "truth"),
                linewidth = 1.5) +
      geom_smooth(aes(colour = "average estimate"),
                  method = "gam", formula = y ~ s(x),
                  se = FALSE, linewidth = 1) +
      geom_line(data = cov_line,
                aes(x = .data[[time_axis]], y = cov_y),
                colour = "blue", linetype = "dotted") +
      geom_hline(yintercept = ref_y, linetype = "dashed", color = "red") +
      facet_wrap(~ slice_lab, ncol = ncol_facets) +
      scale_y_continuous(
        name = y_lab,
        sec.axis = sec_axis(~ (. - rng_min) / rng_span,
                            name   = "coverage",
                            breaks = seq(0, 1, 0.25))
      ) +
      scale_colour_brewer("", palette = "Dark2") +
      labs(x = time_axis,
           title = paste("transition", trans)) +
      theme_bw(base_size = font_size) +
      theme(
        legend.position.inside = c(0.03, 0.97),
        legend.justification   = c("left", "top"),
        legend.background      = element_rect(fill = scales::alpha("white", 0.8),
                                              colour = NA)
      )
  }

  ## ---- return list of plots (one per grouping) ------------------------
  plot_data %>%
    split(.$grouping) %>%
    map(make_plot) %>%
    set_names(unique(plot_data$grouping))
}


# plot_one_fixedEffect <- function(var_name) {

#   df <- res_fe %>%
#     filter(variable == var_name) %>%
#     mutate(
#       model = sub("_[^_]+$", "", model),
#       algo    = factor(model, levels = algo_levels),
#       bs      = factor(bs,        levels = bs_levels),
#       problem = recode(problem, !!!problem_labs)
#     )

#   true_val <- unique(df$true_value)

#   ggplot(df,
#          aes(x = algo,
#              y = estimate,
#              fill = bs,
#              group = interaction(algo, bs))) +

#     geom_boxplot(outlier.shape = NA,
#                  position = position_dodge(width = 0.75)) +

#     geom_hline(yintercept = true_val,
#                colour = "orange",
#                linewidth = 1.3) +

#     scale_fill_brewer(palette = "Dark2", name = "bs") +

#     facet_wrap(~ problem, nrow = 1) +

#     labs(
#       x     = "algorithm",
#       y     = var_name,
#       title = paste("Transition:", var_name)
#     ) +

#     theme_bw(base_size = 14) +
#     theme(
#       axis.text.x    = element_text(angle = 45, hjust = 1),
#       panel.spacing.x = unit(1.2, "lines")
#     )
# }

plot_fixedEffects_boxplots <- function(res_fe, algo_levels = c("algo_timeScales","algo_stratified"),
                                      bs_levels = c("ps","fs"), font_size = 14) {

  # Helper function to create facet labels as plotmath strings
  parse_facet_label <- function(var_string) {
    # Splits "beta_1_01" into "1" and "01"
    parts <- str_split(var_string, "_", simplify = TRUE)
    x_sub <- parts[2]
    last_part <- parts[3]
    # Splits "01" into "0" and "1"
    other_subs <- paste(str_split(last_part, "", simplify = TRUE), collapse = ",")

    # Creates a string that label_parsed can understand.
    # The list() function groups the subscripts correctly.
    paste0("beta[list(x[", x_sub, "],", other_subs, ")]")
  }

  df <- res_fe %>%
    mutate(
      problem_label = recode(problem,
                             "sim_stratified_fe" = "SSTS DGP",
                             "sim_timeScales_fe" = "MTS DGP"),
      model = sub("_[^_]+$", "", model),
      algo  = factor(model, levels = algo_levels),
      bs    = factor(bs,    levels = bs_levels),
      bs_label = recode(bs, ps = "Penalized Splines", fs = "Factor Smooth"),
      algo_label = recode(as.character(algo),
                          "algo_timeScales" = "MTS PAM",
                          "algo_stratified" = "SSTS PAM"),

      x_label = factor(
        paste0(problem_label, " &\n", algo_label),
        levels = c(
          "SSTS DGP &\nSSTS PAM",
          "SSTS DGP &\nMTS PAM",
          "MTS DGP &\nSSTS PAM",
          "MTS DGP &\nMTS PAM"
        )
      ),

      # --- THIS IS THE CORRECTED FACET LABEL CREATION ---
      # Use map_chr to apply the function to each element of 'variable'
      facet_label = purrr::map_chr(variable, parse_facet_label)
    )

  # truth_df needs the new facet_label to map correctly
  truth_df <- df %>%
    distinct(facet_label, true_value) %>%
    rename(true_val = true_value)

  # Ensure the facet labels in the plot are ordered correctly based on the original variable
  df$facet_label <- factor(df$facet_label, levels = unique(df$facet_label[order(df$variable)]))

  ggplot(df,
         aes(x = x_label, y = estimate,
             fill = bs_label,
             color = algo_label
             )) +
    geom_boxplot(outlier.shape = NA,
                 position = position_dodge(width = 0.85),
                 linewidth = 0.6) +
    geom_hline(data = truth_df,
               aes(yintercept = true_val),
               colour = "darkorange", linewidth = 1.1) +

    scale_fill_manual(
      name = "Smooth Type",
      values = c("Penalized Splines" = "#66C2A5", "Factor Smooth" = "#FC8D62")
    ) +
    scale_color_manual(
      name = "Algorithm",
      values = c("SSTS PAM" = "#A6A6A6", "MTS PAM" = "#4F4F4F"),
      breaks = c("SSTS PAM", "MTS PAM")
    ) +

    # --- USE THE NEW FACET LABEL WITH label_parsed ---
    facet_wrap(~ facet_label, ncol = 2, labeller = label_parsed) +
    labs(x = NULL, y = "Estimate") +
    guides(
      fill = guide_legend(title.position="top", title.hjust = 0.5),
      color = guide_legend(title.position="top", title.hjust = 0.5)
    ) +
    theme_bw(base_size = font_size) +
    theme(
      axis.text.x      = element_text(hjust = 0.5),
      panel.spacing    = unit(1.2, "lines"),
      plot.title       = element_blank(),
      strip.background = element_rect(fill = "gray92", color = "gray92"),
      strip.text       = element_text(size = 14, margin = margin(t = 5, b = 5)),
      legend.position  = "bottom",
      legend.box       = "horizontal",
      legend.title     = element_text(size = 11, face = "bold"),
      legend.margin    = margin(t = 10, b = 0)
    )
}


plot_fixedEffects_coverages <- function(
  res_fe,
  algo_levels = c("algo_timeScales","algo_stratified"),
  bs_levels   = c("ps","fs")
) {

  df <- res_fe %>%
    mutate(
      model = sub("_[^_]+$", "", model),
      algo  = factor(model, levels = algo_levels),
      bs    = factor(bs,    levels = bs_levels)
    )

  # summarise coverage + exact binomial 95% CI
  cov_sum <- df %>%
    group_by(variable, algo, bs) %>%
    summarise(
      n        = dplyr::n(),
      x        = sum(coverage, na.rm = TRUE),
      coverage = x / n,
      .groups  = "drop"
    ) %>%
    mutate(
      ci   = map2(x, n, ~ stats::binom.test(.x, .y)$conf.int),
      lower = map_dbl(ci, 1),
      upper = map_dbl(ci, 2)
    ) %>%
    select(-ci)

  dodge <- position_dodge(width = 0.7)

  ggplot(cov_sum,
         aes(x = algo, y = coverage,
             colour = bs, group = interaction(algo, bs))) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  position = dodge, width = 0.15) +
    geom_point(position = dodge, size = 1.5) +
    geom_hline(yintercept = 0.95, colour = "red",
               linetype = "dashed", linewidth = 1) +
    scale_colour_brewer(palette = "Dark2", name = "bs") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = scales::percent_format(accuracy = 1)
    ) +
    facet_wrap(~ variable, ncol = 2) +   # 2×2 grid over the 4 variables
    labs(x = NULL, y = "coverage") +
    theme_bw(base_size = 14) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      panel.spacing = unit(1.2, "lines"),
      plot.title = element_blank(),
      strip.text = element_text(size = 11)
    )
}


plot_one_bh_coverage <- function(df, transition_name) {

  problem_labs   <- c(sim_timeScales_bh = "time-scales DGP",
                      sim_stratified_bh = "stratified DGP")
  metric_levels  <- c("log-hazard", "cumu hazard", "trans prob")

  df_base <- df %>%
    dplyr::filter(transition == transition_name) %>%
    dplyr::mutate(
      algo = sub("_[^_]+$", "", model),
      bs   = sub(".*_", "", model)
    )

  df_cov <- df_base %>%
    tidyr::pivot_longer(
      cols = c(avg_loghazard_coverage,
               avg_cumu_hazard_coverage,
               avg_trans_prob_coverage),
      names_to  = "metric",
      values_to = "coverage"
    ) %>%
    dplyr::mutate(
      metric = dplyr::recode(metric,
        "avg_loghazard_coverage"   = "log-hazard",
        "avg_cumu_hazard_coverage" = "cumu hazard",
        "avg_trans_prob_coverage"  = "trans prob"),
      metric  = factor(metric, levels = metric_levels),
      algo    = factor(algo, levels = c("algo_timeScales", "algo_stratified")),
      bs      = factor(bs,   levels = c("ps", "fs")),
      problem = dplyr::recode(problem, !!!problem_labs)
    )

  ci_lower_cols <- c("avg_loghazard_coverage_lower",
                     "avg_cumu_hazard_coverage_lower",
                     "avg_trans_prob_coverage_lower")
  ci_upper_cols <- c("avg_loghazard_coverage_upper",
                     "avg_cumu_hazard_coverage_upper",
                     "avg_trans_prob_coverage_upper")
  have_ci <- all(ci_lower_cols %in% names(df)) && all(ci_upper_cols %in% names(df))

  if (have_ci) {
    df_low <- df_base %>%
      tidyr::pivot_longer(dplyr::all_of(ci_lower_cols),
                          names_to = "metric", values_to = "lower") %>%
      dplyr::mutate(
        problem = dplyr::recode(problem, !!!problem_labs),
        metric  = dplyr::recode(metric,
          "avg_loghazard_coverage_lower"   = "log-hazard",
          "avg_cumu_hazard_coverage_lower" = "cumu hazard",
          "avg_trans_prob_coverage_lower"  = "trans prob"),
        metric = factor(metric, levels = metric_levels)
      )

    df_up <- df_base %>%
      tidyr::pivot_longer(dplyr::all_of(ci_upper_cols),
                          names_to = "metric", values_to = "upper") %>%
      dplyr::mutate(
        problem = dplyr::recode(problem, !!!problem_labs),
        metric  = dplyr::recode(metric,
          "avg_loghazard_coverage_upper"   = "log-hazard",
          "avg_cumu_hazard_coverage_upper" = "cumu hazard",
          "avg_trans_prob_coverage_upper"  = "trans prob"),
        metric = factor(metric, levels = metric_levels)
      )

    df_cov <- df_cov %>%
      dplyr::left_join(df_low, by = c("transition","problem","model","algo","bs","metric")) %>%
      dplyr::left_join(df_up,  by = c("transition","problem","model","algo","bs","metric")) %>%
      dplyr::mutate(metric = factor(metric, levels = metric_levels))   # <- enforce order again
  }

  dodge <- position_dodge(width = 0.7)

  ggplot(df_cov, aes(x = algo, y = coverage, colour = bs, group = interaction(algo, bs))) +
    { if (have_ci) geom_errorbar(aes(ymin = lower, ymax = upper), position = dodge, width = 0.15) else NULL } +
    geom_point(position = dodge, size = 1.5) +
    facet_grid(problem ~ metric, switch = "y") +
    scale_colour_brewer(palette = "Dark2", name = "bs") +
    geom_hline(yintercept = 0.95, colour = "red", linetype = "dashed", linewidth = 1) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "average coverage",
         title = paste("Transition", transition_name)) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.spacing = unit(1, "lines"),
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0))
}
create_bias_rmse_plot <- function(df_lines, m, title, df_hist, scale_f, font_size = 14) {

  # --- Data Preparation ---
  df_hist_labeled <- df_hist %>%
    mutate(
      DGP = recode(problem,
                   "sim_timeScales_bh" = "MTS DGP",
                   "sim_stratified_bh" = "SSTS DGP")
    )

  df_lines_labeled <- df_lines %>%
    mutate(
      `DGP-by-model combination` = recode(
        problem_model,
        "sim_stratified_bh / algo_stratified_ps" = "SSTS DGP & STSS PAM (PS)",
        "sim_stratified_bh / algo_stratified_fs" = "SSTS DGP & STSS PAM (FS)",
        "sim_timeScales_bh / algo_timeScales_ps" = "MTS DGP & MTS PAM (PS)",
        "sim_timeScales_bh / algo_timeScales_fs" = "MTS DGP & MTS PAM (FS)",
        "sim_stratified_bh / algo_timeScales_ps" = "SSTS DGP & MTS PAM (PS)",
        "sim_stratified_bh / algo_timeScales_fs" = "SSTS DGP & MTS PAM (FS)",
        "sim_timeScales_bh / algo_stratified_ps" = "MTS DGP & SSTS PAM (PS)",
        "sim_timeScales_bh / algo_stratified_fs" = "MTS DGP & SSTS PAM (FS)"
      )
      # Factor ordering can be added here if needed
    )

  # --- Plotting ---
  ggplot() +
    ## Histogram
    geom_col(data = df_hist_labeled,
             aes(x = tend, y = value_scaled, fill = DGP),
             position = "identity", alpha = 0.35, width = 0.95) +

    ## Bias / RMSE lines
    geom_line(data = df_lines_labeled %>% filter(metric == m),
              aes(x = tend, y = value,
                  colour = `DGP-by-model combination`,
                  group = `DGP-by-model combination`)) +
    geom_point(data = df_lines_labeled %>% filter(metric == m),
               aes(x = tend, y = value,
                   colour = `DGP-by-model combination`,
                   group = `DGP-by-model combination`),
               size = 1) +

    # Facets and Labels
    facet_wrap(~ transition, nrow = 1,
               labeller = labeller(transition = c("1" = "Transition 1", "2" = "Transition 2", "3" = "Transition 3"))) +
    labs(x = "Time", y = NULL, title = title,
         fill = "DGP", colour = "DGP-by-model \ncombination") +

    # Axes and Scales
    scale_y_continuous(
      sec.axis = sec_axis(~ . / scale_f, name = "Relative Event Frequency (in %)")
    ) +
    scale_fill_manual(values = c("MTS DGP" = "#1f77b4", "SSTS DGP" = "#e377c2")) +

    # Guides and Theme
    guides(
      fill   = guide_legend(order = 1, override.aes = list(alpha = 0.35)),
      colour = guide_legend(order = 2, nrow = 4, ncol = 2, byrow = TRUE)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      axis.title      = element_text(size = font_size),
      axis.text       = element_text(size = font_size),
      # --- MODIFICATION HERE ---
      # Add margin to legend titles. `b` is bottom, `r` is right.
      legend.title    = element_text(size = font_size, margin = margin(b = 5, r = 15)),
      legend.text     = element_text(size = font_size),
      strip.text      = element_text(size = font_size),
      plot.title      = element_text(size = font_size, hjust = 0.5)
    )
}


generate_coverage_latex_bh <- function(summary_df) {

  escape_latex <- function(text) {
    if (is.na(text)) return("")
    text %>%
      str_replace_all("&", "\\\\&") %>% str_replace_all("%", "\\\\%") %>%
      str_replace_all("_", "\\\\_")
  }

  data_cols <- names(summary_df)[-c(1, 2)]
  col_info <- tibble(
    dgp = sapply(strsplit(data_cols, "_"), `[`, 1),
    model = sapply(strsplit(data_cols, "_"), `[`, 2),
    smoother = sapply(strsplit(data_cols, "_"), `[`, 3)
  )

  dgp_levels <- c("SSTS DGP", "MTS DGP")
  model_levels <- c("SSTS PAM", "MTS PAM")

  dgp_structure <- col_info %>%
    mutate(dgp = factor(dgp, levels = dgp_levels)) %>%
    group_by(dgp) %>% summarise(n = n(), .groups = 'drop') %>% arrange(dgp)

  model_structure <- col_info %>%
    mutate(dgp = factor(dgp, levels = dgp_levels), model = factor(model, levels = model_levels)) %>%
    group_by(dgp, model) %>% summarise(n = n(), .groups = 'drop') %>% arrange(dgp, model)

  # --- CHANGE 1: Abbreviate the DGP headers ---
  dgp_headers <- mapply(function(dgp, n) { paste0("\\multicolumn{", n, "}{c}{", dgp, "}") }, dgp_structure$dgp, dgp_structure$n)
  dgp_header_line <- paste("& &", paste(dgp_headers, collapse = " & "), "\\\\")

  model_headers <- mapply(function(model, n) { paste0("\\multicolumn{", n, "}{c}{", model, "}") }, model_structure$model, model_structure$n)
  model_header_line <- paste("& &", paste(model_headers, collapse = " & "), "\\\\")

  smoother_headers_multiline <- str_replace_all(col_info$smoother,
    c("Penalized Splines" = "\\\\shortstack[c]{Penalized\\\\\\\\Splines}",
      "Factor Smooth" = "\\\\shortstack[c]{Factor\\\\\\\\Smooth}")
  )
  smoother_header_line <- paste(
    "Transition & \\shortstack[l]{Quantity of \\\\ Interest}",
    paste(smoother_headers_multiline, collapse = " & "),
    sep = " & "
  )

  cmid_start_col <- 3
  dgp_cmidrules <- c()
  for (n in dgp_structure$n) {
    cmid_end_col <- cmid_start_col + n - 1
    dgp_cmidrules <- c(dgp_cmidrules, paste0("\\cmidrule(lr){", cmid_start_col, "-", cmid_end_col, "}"))
    cmid_start_col <- cmid_end_col + 1
  }

  cmid_start_col <- 3
  model_cmidrules <- c()
  for (n in model_structure$n) {
    cmid_end_col <- cmid_start_col + n - 1
    model_cmidrules <- c(model_cmidrules, paste0("\\cmidrule(lr){", cmid_start_col, "-", cmid_end_col, "}"))
    cmid_start_col <- cmid_end_col + 1
  }

  num_total_data_cols <- length(data_cols)
  tabular_def <- paste0("\\begin{tabular}{lc", paste(rep("c", num_total_data_cols), collapse = ""), "}")

  header <- c(
    "\\begin{table*}[htbp]", "\\centering", "\\caption{\\captiontsbhcoverage}",
    "\\label{tab:sim-ts-bh-coverage}", "\\begin{sideways}",
    # --- CHANGE 2: Further reduce inter-column space ---
    "\\setlength{\\tabcolsep}{2.5pt}",
    tabular_def, "\\toprule",
    dgp_header_line, paste(dgp_cmidrules, collapse = " "),
    model_header_line, paste(model_cmidrules, collapse = " "),
    paste0(smoother_header_line, " \\\\"),
    "\\midrule"
  )

  body_rows <- c()
  for (i in 1:nrow(summary_df)) {
    transition_label <- escape_latex(summary_df$Transition[i]) %>% str_replace_all("->", " $\\\\rightarrow$ ")
    quantity_label <- paste0("\\makecell[l]{", escape_latex(summary_df$`Quantity of Interest`[i]), "}") %>% str_replace_all("\\\\n", "\\\\\\\\")

    row_labels <- c(transition_label, quantity_label)

    data_cells <- summary_df[i, -c(1, 2)]
    formatted_cells <- sapply(data_cells, function(cell_content) {
      if (is.na(cell_content)) return("")
      else {
        content_with_break <- str_replace(cell_content, " \\(", "\\\\\\\\(")
        return(paste0("\\makecell[t]{", content_with_break, "}"))
      }
    })

    full_row <- paste(c(row_labels, formatted_cells), collapse = " & ")
    body_rows <- c(body_rows, paste0(full_row, " \\\\"))
  }

  footer <- c("\\bottomrule", "\\end{tabular}", "\\end{sideways}", "\\end{table*}")
  full_latex_code <- paste(c(header, body_rows, footer), collapse = "\n")
  return(full_latex_code)
}


generate_coverage_latex_fe <- function(summary_df, rounding_digits = 2) {

  # --- Step 1: Internal Data Preparation and Pivoting ---
  summary_formatted <- summary_df %>%
    mutate(
      dgp_label = recode(problem,
                         "sim_stratified_fe" = "SSTS DGP",
                         "sim_timeScales_fe" = "MTS DGP"),
      model_label = recode(model_name,
                           "algo_stratified" = "SSTS PAM",
                           "algo_timeScales" = "MTS PAM"),
      smoother_label = recode(bs_type,
                              "ps" = "Penalized Splines",
                              "fs" = "Factor Smooth"),
      cell_value = sprintf(paste0("%.", rounding_digits, "f (%.", rounding_digits, "f; %.", rounding_digits, "f)"),
                           coverage, coverage_lower, coverage_upper)
    )

  # Define the order for the final columns
  dgp_order <- c("SSTS DGP", "MTS DGP")
  model_order <- c("SSTS PAM", "MTS PAM")
  smoother_order <- c("Penalized Splines", "Factor Smooth")

  # Create all possible column names in the correct final order
  final_col_order <- expand_grid(dgp = dgp_order, model = model_order, smoother = smoother_order) %>%
    unite("col_name", dgp, model, smoother, sep = "_") %>%
    pull(col_name)

  wide_df <- summary_formatted %>%
    pivot_wider(
      id_cols = variable,
      names_from = c(dgp_label, model_label, smoother_label),
      names_sep = "_",
      values_from = cell_value
    ) %>%
    rename("Transition" = variable) %>%
    mutate(Transition = str_replace_all(Transition,
      c("beta_1_01" = "$\\\\beta_{x_1,0,1}$",
        "beta_1_03" = "$\\\\beta_{x_1,0,3}$",
        "beta_1_12" = "$\\\\beta_{x_1,1,2}$",
        "beta_1_13" = "$\\\\beta_{x_1,1,3}$")
    )) %>%
    # Use the ordered vector to select and arrange columns
    select(any_of(c("Transition", final_col_order)))


  # --- Step 2: LaTeX Generation ---

  escape_latex <- function(text) {
    if (is.na(text)) return("")
    return(text)
  }

  data_cols <- names(wide_df)[-1]
  col_info <- tibble(
    dgp = sapply(strsplit(data_cols, "_"), `[`, 1),
    model = sapply(strsplit(data_cols, "_"), `[`, 2),
    smoother = sapply(strsplit(data_cols, "_"), `[`, 3)
  )

  dgp_levels <- c("SSTS DGP", "MTS DGP")
  model_levels <- c("SSTS PAM", "MTS PAM")

  dgp_structure <- col_info %>%
    mutate(dgp = factor(dgp, levels = dgp_levels)) %>%
    group_by(dgp) %>% summarise(n = n(), .groups = 'drop') %>% arrange(dgp)

  model_structure <- col_info %>%
    mutate(dgp = factor(dgp, levels = dgp_levels), model = factor(model, levels = model_levels)) %>%
    group_by(dgp, model) %>% summarise(n = n(), .groups = 'drop') %>% arrange(dgp, model)

  dgp_headers <- mapply(function(dgp, n) { paste0("\\multicolumn{", n, "}{c}{", dgp, "}") }, dgp_structure$dgp, dgp_structure$n)
  dgp_header_line <- paste("&", paste(dgp_headers, collapse = " & "), "\\\\")

  model_headers <- mapply(function(model, n) { paste0("\\multicolumn{", n, "}{c}{", model, "}") }, model_structure$model, model_structure$n)
  model_header_line <- paste("&", paste(model_headers, collapse = " & "), "\\\\")

  smoother_headers_multiline <- str_replace_all(col_info$smoother,
    c("Penalized Splines" = "\\\\shortstack[c]{Penalized\\\\\\\\Splines}",
      "Factor Smooth" = "\\\\shortstack[c]{Factor\\\\\\\\Smooth}")
  )
  smoother_header_line <- paste("Transition", paste(smoother_headers_multiline, collapse = " & "), sep = " & ")

  cmid_start_col <- 2
  dgp_cmidrules <- c()
  for (n in dgp_structure$n) {
    cmid_end_col <- cmid_start_col + n - 1
    dgp_cmidrules <- c(dgp_cmidrules, paste0("\\cmidrule(lr){", cmid_start_col, "-", cmid_end_col, "}"))
    cmid_start_col <- cmid_end_col + 1
  }

  cmid_start_col <- 2
  model_cmidrules <- c()
  for (n in model_structure$n) {
    cmid_end_col <- cmid_start_col + n - 1
    model_cmidrules <- c(model_cmidrules, paste0("\\cmidrule(lr){", cmid_start_col, "-", cmid_end_col, "}"))
    cmid_start_col <- cmid_end_col + 1
  }

  num_total_data_cols <- length(data_cols)
  tabular_def <- paste0("\\begin{tabular}{l", paste(rep("c", num_total_data_cols), collapse = ""), "}")

  header <- c(
    "\\begin{table*}[htbp]", "\\centering", "\\caption{\\captionfecoverage}",
    "\\label{tab:sim-ts-fe-coverage}",
    "\\begin{sideways}",
    "\\setlength{\\tabcolsep}{3pt}",
    tabular_def, "\\toprule",
    dgp_header_line, paste(dgp_cmidrules, collapse = " "),
    model_header_line, paste(model_cmidrules, collapse = " "),
    paste0(smoother_header_line, " \\\\"),
    "\\midrule"
  )

  body_rows <- c()
  for (i in 1:nrow(wide_df)) {
    transition_label <- escape_latex(wide_df$Transition[i])

    data_cells <- wide_df[i, -1]
    formatted_cells <- sapply(data_cells, function(cell_content) {
      if (is.na(cell_content)) return("")
      else {
        content_with_break <- str_replace(cell_content, " \\(", "\\\\\\\\(")
        return(paste0("\\makecell[t]{", content_with_break, "}"))
      }
    })

    full_row <- paste(c(transition_label, formatted_cells), collapse = " & ")
    body_rows <- c(body_rows, paste0(full_row, " \\\\"))
  }

  footer <- c("\\bottomrule", "\\end{tabular}", "\\end{sideways}", "\\end{table*}")
  full_latex_code <- paste(c(header, body_rows, footer), collapse = "\n")
  return(full_latex_code)
}


