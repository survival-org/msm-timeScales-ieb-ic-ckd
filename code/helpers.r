library(data.table)
library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)

read_file <- function(a_infile, header="auto",sep="auto",fill=FALSE, ...){
	tblInRaw <- fread(a_infile, header=header, sep=sep)
	tblIn <- as.data.frame(tblInRaw,fill=fill)
	return(tblIn)
}

sim_wrapper <- function(data, job, formula, n = 500, time_grid = seq(0, 10, by = 0.05)) {

  # create data set with covariates
  df <- tibble::tibble(x1 = runif(n, -3, 3), x2 = runif(n, 0, 6))

  ndf <- sim_pexp(
    formula = formula,
    data = df, cut = time_grid)

  ndf <- add_interval_censoring(ndf)

  return(ndf)
}

ci_wrapper <- function(
  data,
  job,
  instance,
  bs      = "ps",
  k       = 10,
  ci_type = "default") {

  # instance <- sim_wrapper()
  ped <- as_ped(
    data    = instance,
    formula = Surv(time, status) ~ x1 + x2,
    id      = "id")

  form <- paste0("ped_status ~ s(tend, bs='", bs, "', k=", k, ") + s(x1) + s(x2)")

  mod <- gam(
    formula = as.formula(form),
    data    = ped,
    family  = poisson(),
    offset  = offset,
    method  = "REML")

  f0 <- function(t) {dgamma(t, 8, 2) * 6}

  # create new data set
  nd <- make_newdata(ped, tend = unique(tend), x1 = c(0), x2 = c(3)) %>% # tbd: check x1 and x2 values here! (currently: mean values)
    mutate(
      true_hazard = exp(-3.5 + f0(tend) -0.5 * x1 + sqrt(x2)),
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
    summarize_all(mean) %>%
    mutate(method = "direct")

    return(nd)

}

add_interval_censoring <- function()