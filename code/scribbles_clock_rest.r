library(dplyr)
library(flexsurv)
library(eha)
library(survival)

data("infants", package = "eha")
infants <- infants %>% mutate(time = exit - enter)
head(infants)


cox_leftTrunc <- coxph(
  formula = Surv(enter, exit, event) ~ mother + age + sex + parish + civst + ses,
  data    = infants)
summary(cox_leftTrunc)

cox_time <- coxph(
    formula = Surv(time, event) ~ mother + age + sex + parish + civst + ses,
    data    = infants
)
summary(cox_time)

cox_time_adj <- coxph(
    formula = Surv(time, event) ~ mother + age + sex + parish + civst + ses + enter,
    data    = infants
)
summary(cox_time_adj)

cox_time_adj_smooth <- coxph(
    formula = Surv(time, event) ~ mother + age + sex + parish + civst + ses + pspline(enter),
    data    = infants
)
summary(cox_time_adj_smooth)

cbind(
        cox_time_adj      = exp(coef(cox_time_adj))[1:7],
    cox_time_adj_smooth = c(exp(coef(cox_time_adj_smooth))[1:6],NA),
  cox_leftTrunc     = c(exp(coef(cox_leftTrunc))[1:6],NA),
  cox_time          = c(exp(coef(cox_time))[1:6],NA))

