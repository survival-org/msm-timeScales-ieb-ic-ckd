
# covariates (separate effects being estimated for each transition):
# - age at first timepoint in a given state
# - age at last timepoint in a given state (or alternatively, time since first timepoint in that state)
# - age at measurement timepoints in-between?
# - eGFR at first timepoint in a given state
# - eGFR at last timepoint in a given state, further eGFR measurements? (could mask genetic effects)
# - sex
# - diabetes status
# - genetics: which one? 12+11? 13+11?

# non-linearities, time-varying effects:
# - age: non-linear, TVE? -> 2-dimensional tensor splines?
# - eGFR: non-linear, TVE? -> 2-dimensional tensor splines?
# - sex: TVE?
# - diabetes: TVE?
# - genetics: linear (or 0-1-2?), TVE?, interactions?

# load packages ----
library(dplyr)
library(data.table)
library(etm)
library(mvna)
library(kmi)
# devtools::load_all("/wis37138/kidney_survival/code/pammtools")
devtools::load_all("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/pammtools") # before doing this: git fetch origin multi-state; git checkout multi-state

# set directories ----
# dir_data <- "/wis37138/kidney_survival/data/ukb"
# dir_out <- "/wis37138/kidney_survival/results/ukb"
# file_events_ukb <- "events_ukb_noAKI6months_ni10+.txt"
setwd("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction")
dir_data <- "data/ukb"
dir_out <- "results/ukb"
file_pheno_ukb <- "pheno_ukb_noAKI6months_ni10+.txt"
file_events_ukb <- "events_ukb_noAKI6months_ni10+.txt"

# load data ----
pheno_ukb <- fread(file.path(dir_data, file_pheno_ukb))
events_ukb <- fread(file.path(dir_data, file_events_ukb))

# Aalen-Johansen estimator ----
events_ukb_aj <- events_ukb %>%
    rename(entry = tstart, exit = tstop2) %>%
    select(id, from, to, entry, exit)

## build transition matrix
tra.idm <- matrix(FALSE, 4, 4, dimnames = list(c(0, 1, 2, 3), c(0, 1, 2, 3)))
tra.idm[1, 2:4] <- TRUE
tra.idm[2, 3:4] <- TRUE
tra.idm[3, 4] <- TRUE

## mvna
mvna.idm <- mvna(events_ukb_aj, c("0", "1", "2", "3"), tra.idm, "cens")

if (require(lattice)){
  xyplot(mvna.idm, xlim=c(0, 100))
}

plot(mvna.idm$'0 1'$time, mvna.idm$'0 1'$na, xlim=c(0,100))

## etm
### fit
etm.idm <- etm(events_ukb_aj, c("0", "1", "2", "3"), tra.idm, "cens", s = 0)

### plot
transitions <- c("0 1", "0 2", "1 2", "1 3", "2 3")
  par(mfrow=c(2,3))
  for(i in 1:length(transitions)){
    plot(etm.idm, tr.choice = transitions[i], conf.int = TRUE,
        lwd = 2, legend = FALSE, ylim = c(0, 1),
        xlim = c(0, 100), xlab = "Years",
        ci.fun = "cloglog", main = transitions[i])
  }
  par(mfrow=c(1,1))

### table
# summary(etm.idm)$"0 1"[, c("P", "lower", "upper")]

### landmark method for probability plots (given from status at time s, probability of being in state to at time t, s <= t)
time.points <- TBD
landmark.etm <- lapply(time.points, function(start.time) {
    etm(events_ukb_aj, c("0", "1", "2"), tra.idm, "cens", start.time)
    })

#### plot landmark
TBD

# PEM & PAM without covariates ----
events_ukb_ct <- events_ukb %>% add_counterfactual_transitions()

ped_events_ukb <- as_ped_multistate(
  data       = events_ukb_ct,
  formula    = Surv(tstart, tstop2, status) ~ 1,
  transition = "transition",
  id         = "id",
  censor_code = 0,
  timescale  = "calendar")

dim(events_ukb)
dim(events_ukb_ct)
dim(ped_events_ukb)

# predict hazards
# as factor: transition in cal_icu.pneu
# ggf in as_ped gleich als faktorvariable ausgeben.
pam <- mgcv::gam(
  formula = ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition),
  data = ped_events_ukb,
  family=poisson(),
  offset=offset)
summary(pam)

plot(pam, xlim = c(0,50), ylim = c(-20, 20), page=1)
plot.gam(pam, select = 4, ylim=c(-1,1))

meeting_test <- make_newdata(cal_icu.pneu, tend = unique(tend), transition=unique(transition), age = quantile(age, probs=c(0.05, 0.95))) %>%
  group_by(transition, age) %>%
  add_cumu_hazard(pam)

# use add_trans_prob to calculate transition prob
test <- meeting_test |> add_trans_prob(pam)
ggplot(test, aes(x=tend, y=trans_prob)) +
  geom_line(aes(col=as.factor(age))) +
  facet_wrap(~transition, ncol = 1) +
  xlim(c(0, 100)) +
  ylim(c(0,1))
