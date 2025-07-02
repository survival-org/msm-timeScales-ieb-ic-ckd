library(survival)
library(dplyr)
devtools::load_all("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/pammtools")

# example from pammtools ----
# set number of observations/subjects
n <- 250
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
 as_tibble()
# the formula which specifies how covariates affet the hazard rate
f0 <- function(t) {
 dgamma(t, 8, 2) *6
}
plot(x=1:10, y=f0(1:10), type="l")
form <- ~ -3.5 + f0(t) -0.5*x1 + sqrt(x2)
set.seed(24032018)
sim_df <- sim_pexp(form, df,  1:10)
head(sim_df)
plot(survfit(Surv(time, status)~1, data = sim_df ))

# for control, estimate with Cox PH
mod <- coxph(Surv(time, status) ~ x1 + pspline(x2), data=sim_df)
coef(mod)[1]
layout(matrix(1:2, nrow=1))
termplot(mod, se = TRUE)

# and using PAMs
layout(1)
ped <- sim_df %>% as_ped(Surv(time, status)~., max_time=10)
library(mgcv)
pam <- gam(ped_status ~ s(tend) + x1 + s(x2), data=ped, family=poisson, offset=offset)
coef(pam)[2]
plot(pam, page=1)

# CR (Giovanni) ----

# define true functions for plotting
f0 <- function(t) {
  dgamma(t, 8, 2) * 12 - 0.5
}

fx <- function(x, transition) {
  (transition == "1->2")*(0.5 + 0.25*x**3) + (transition == "2->3")*(0.4 - x**2)
}

ft <- function(t, transition) {
  (transition == "1->3")*f0(t)
}

# # define hazard for transitions 1->2 & 1->3
# form_stepone <- formula( ~ 0.5 + 0.05*x1**3 + f0(t) | 0.6 - 0.05*x1**2  )
# # define hazard for transition 2->3
# form_steptwo <- formula( ~ - 0.4 + sin(x2) + (t+x3 <= 3) * f0(t + x3))

# define hazard for transitions 1->2 & 1->3
form_stepone <- formula( ~ 0.5 + 0.25*x1**3| f0(t)  )
# define hazard for transition 2->3
form_steptwo <- formula( ~ 0.4 - x1**2)

nSim = 20

n = 2500
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  from = 1,
  t = 0)

cut =  seq(0, 3, by = 0.01)
seq_x = seq(-2, 2, by = 0.1)
seq_t = seq(0, 3, by = 0.1)

i = 1

for (i in 1:nSim) {
  sim_df_stepone <- pammtools:::sim_pexp_cr(form_stepone, data, cut) %>%
    mutate(from = 1
           , to = type + 1
           , transition = case_when(
             type == 1 ~ "1->2",
             type == 2 ~ "1->3",
             .default = "err"
           )) %>%
    rename(tstart = t, tstop = time) %>%
    filter(status == 1)

  # second step
  # take all relevant subjects, i.e. subjects with transition 1->2, i.e. currently in state 2
  data_steptwo <- sim_df_stepone %>%
    filter(to == 2, status == 1) %>%
    select(id, x1, from, t = tstop)

  sim_df_steptwo <- pammtools:::sim_pexp_cr(form_steptwo, data_steptwo, cut) %>%
    mutate(from = 2
           , to = type + 2
           , transition = case_when(
             type == 1 ~ "2->3",
             .default = "err"
           )
        #    , status = ifelse(time + t > 3, 0, status)
        #    , time = min(time + t, 3)
           , hazard2 = 0) %>%
    rename(tstart = t, tstop = time) %>%
    filter(status == 1)

  sim_df <- rbind(sim_df_stepone, sim_df_steptwo)
  sim_df <- sim_df %>% add_counterfactual_transitions()
  #
  # test <- sim_df %>% mutate(time = ifelse(tstop > 2, 1, 0)) %>% filter(status == 1)
  # table(test$time, test$transition)

  # go on with pamm tools procedure
  cal_sim_df <- as_ped_multistate(
    data       = sim_df,
    formula    = Surv(tstart, tstop, status)~ transition + x1,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")
  #
  # dim(cal_sim_df)
  # head(cal_sim_df)

  ctrl <- gam.control(trace = TRUE)
  bam_sim <- mgcv::bam(ped_status ~ s(tend, by=as.factor(transition))
                       + as.factor(transition)
                       + s(x1, by=as.factor(transition))
                       , data = cal_sim_df
                       , family=poisson()
                       , offset=offset
                       , discrete = T
                       , method = "fREML"
                       , control = ctrl)

  # summary(bam_sim)

  plot(bam_sim, select = 6, xlim = c(-2,2), ylim = c(-1,1))

  if(i == 1) {
    x1_df <- cal_sim_df %>%
      make_newdata(x1 = seq_x, transition = unique(transition)) %>%
      add_term(bam_sim, term = "x1")%>%
      mutate(true.value = fx(x=x1, transition=transition)
             , shift = case_when(
               transition == "1->2" ~ coef(bam_sim)[1]
               , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
               , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
               , TRUE ~ 0
             )
             , fit = fit + shift
             , ci_lower = ci_lower + shift
             , ci_upper = ci_upper + shift
             , iter = i)

    tend_df <- cal_sim_df %>%
      make_newdata(tend = unique(tend), transition = unique(transition)) %>%
      add_term(bam_sim, term = "tend") %>%
      mutate(true.value = ft(tend, transition)
             , shift = case_when(
               transition == "1->2" ~ coef(bam_sim)[1]
               , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
               , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
               , TRUE ~ 0
             )
             , fit = fit + shift
             , ci_lower = ci_lower + shift
             , ci_upper = ci_upper + shift
             , iter = i
      )
  } else {
    x1_df_temp <- cal_sim_df %>%
      make_newdata(x1 = seq_x, transition = unique(transition)) %>%
      add_term(bam_sim, term = "x1") %>%
      mutate(true.value = fx(x=x1, transition=transition)
             , shift = case_when(
               transition == "1->2" ~ coef(bam_sim)[1]
               , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
               , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
               , TRUE ~ 0
             )
             , fit = fit + shift
             , ci_lower = ci_lower + shift
             , ci_upper = ci_upper + shift
             , iter = i)
    x1_df <- rbind(x1_df, x1_df_temp)

    tend_df_temp <- cal_sim_df %>%
      make_newdata(tend = unique(tend), transition = unique(transition)) %>%
      add_term(bam_sim, term = "tend") %>%
      mutate(true.value = ft(tend, transition)
             , shift = case_when(
               transition == "1->2" ~ coef(bam_sim)[1]
               , transition == "1->3" ~ coef(bam_sim)[1] + coef(bam_sim)[2]
               , transition == "2->3" ~ coef(bam_sim)[1] + coef(bam_sim)[3]
               , TRUE ~ 0
             )
             , fit = fit + shift
             , ci_lower = ci_lower + shift
             , ci_upper = ci_upper + shift
             , iter = i
      )
    tend_df <- rbind(tend_df, tend_df_temp)
  }

}

# MSM ----
## setup ----
f0 <- function(t) {
  dgamma(t, 8, 2) * 12 - 0.5
}
plot(x=seq(0,3,0.01), y=f0(seq(0,3,0.01)), type="l")

formulas_list <- list(
    list(0, 1, ~ 0.2 + f0(t)),
    list(0, 4, ~ 0.1 + f0(t) + 0.5*x1**2),
    list(1, 2, ~ 0.1 + 0.25*x1**3 + 0.5*x2 + 0.5*f0(t)),
    list(1, 4, ~ 0.4 + 0.7*x1**2),
    list(2, 3, ~ 0.3 + 0.25*x1**2 + 0.5*x2 + 0.3*x3 + 0.25*f0(t)),
    list(2, 4, ~ f0(t))
)
# formulas_list <- list(
#   list(from = 0, to = 1, formula = ~ -4.5 + f0(t)),
#   list(from = 0, to = "death", formula = ~ -3.0 + f0(t) + 0.5*x1**2),
#   list(1, "death", ~ -2.1 + 0.7*x1**2)
# )

terminal_states <- c(3, 4)
cut =  seq(0, 3, by = 0.01)

n = 2500

set.seed(123)
data <- cbind.data.frame(
  id = 1:n,
  x1 = runif(n, -3, 3),
  x2 = runif(n, 0, 6),
  x3 = runif(n, 0, 3),
  from = 0,
  t = 0)

## only type I right-censoring ----
sim_df_msm <- sim_pexp_msm(formulas_list, data, cut, terminal_states, round = 2, add_counterfactuals = FALSE)
table(sim_df_msm$status)

ped_sim_df_msm <- as_ped_multistate(
    data       = sim_df_msm %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

ctrl <- gam.control(trace = TRUE)
# formula <- ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition) + s(x1, by=as.factor(transition)) +
#                 s(x2, by=as.factor(transition)) + s(x3, by=as.factor(transition))
formula <- ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition) + s(x1, by=as.factor(transition)) + s(x2, by=as.factor(transition))
bam_sim <- mgcv::bam(formula
                    , data = ped_sim_df_msm
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
                    , control = ctrl)
summary(bam_sim)
plot(bam_sim)

## adding type III right-censoring ----
sim_df_msm_rightcens <- sim_df_msm %>%
  add_censoring(type = "right", distribution = "weibull", parameters = c(1, 1), round = 2)

ped_sim_df_msm_rightcens <- as_ped_multistate(
    data       = sim_df_msm_rightcens %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

bam_sim_rightcens <- mgcv::bam(formula
                    , data = ped_sim_df_msm_rightcens
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
                    , control = ctrl)
summary(bam_sim_rightcens)
plot(bam_sim_rightcens)

## adding interval-censoring ----
sim_df_msm_intervalcens <- sim_df_msm %>%
  add_censoring(type = "interval", distribution = "uniform", parameters = NULL, round = 2)

ped_sim_df_msm_intervalcens <- as_ped_multistate(
    data       = sim_df_msm_intervalcens %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

bam_sim_intervalcens <- mgcv::bam(formula
                    , data = ped_sim_df_msm_intervalcens
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
                    , control = ctrl)
summary(bam_sim_intervalcens)
plot(bam_sim_intervalcens)

## adding left-censoring ----
sim_df_msm_leftcens <- sim_df_msm %>%
  add_censoring(type = "left", distribution = "weibull", parameters = c(0.5, 0.5), round = 2)

ped_sim_df_msm_leftcens <- as_ped_multistate(
    data       = sim_df_msm_leftcens %>% add_counterfactual_transitions(),
    formula    = Surv(tstart, tstop, status) ~ .,
    transition = "transition",
    id         = "id",
    censor_code = 0,
    timescale  = "calendar")

bam_sim_leftcens <- mgcv::bam(formula
                    , data = ped_sim_df_msm_leftcens
                    , family=poisson()
                    , offset=offset
                    , discrete = T
                    , method = "fREML"
                    , control = ctrl)
summary(bam_sim_leftcens)
plot(bam_sim_leftcens)
