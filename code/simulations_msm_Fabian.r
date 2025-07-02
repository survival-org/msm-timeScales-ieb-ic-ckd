library(survival)
library(dplyr)
devtools::load_all("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/pammtools")

source("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/msm_kidneyFunction/code/helpers_ukb.r")

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
  t = 20)

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
