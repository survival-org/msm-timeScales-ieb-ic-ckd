library(dplyr)
library(etm)
library(mvna)
library(kmi)
devtools::load_all("C:/Users/ra56yaf/Desktop/Projects/StaBLab/Survival Analysis/survival_kidneyFunction/pammtools") # before doing this: git fetch origin multi-state; git checkout multi-state


# DATAPREP ----

# load data used in BEYERSMANN
data(icu.pneu)

# BEYERSMANN
my.icu.pneu <- icu.pneu %>%
  group_by(id) %>%
  mutate(event = event
         , to2  = ifelse(status == 1, as.character(event), ifelse((n() == 1 & status == 0) | (pneu == 1 & status == 0), "cens", "1"))
         , to   = ifelse((n() == 1 & status == 0) | (pneu == 1 & status == 0), "cens", as.character(ifelse(status == 1, min(event, 2), 1) ))
         , transition = paste0(as.character(pneu), "->", as.character(to))
  ) %>%
  rename(entry = start
         , exit = stop
         , from = pneu
  ) %>%
  select(id
         , entry
         , exit
         , from
         , to
         , to2
         , transition
         , age
         , sex
         # , adm.cens.exit
         # , status
         # , event
  )

# RECALCULATE BEYERSMANN ----

# build transition matrix
tra.idm <- matrix(FALSE, 3, 3, dimnames = list(c(0, 1, 2), c(0, 1, 2)))
tra.idm[1, 2:3] <- TRUE
tra.idm[2, 3] <- TRUE

# print result of transition matrix
print(tra.idm)

mvna.idm <- mvna(my.icu.pneu, c("0", "1", "2"), tra.idm, "cens")

# re-create plots of book
if (require(lattice)){
  xyplot(mvna.idm, xlim=c(0, 100))
}
# same graphs as in Beyersmann --> correct

plot(mvna.idm$'0 1'$time, mvna.idm$'0 1'$na, xlim=c(0,100))

# estimate transition probabilites
etm.idm <- etm(my.icu.pneu, c("0", "1", "2"), tra.idm, "cens", s = 0)

transitions <- c("0 1", "0 2", "1 2")
par(mfrow=c(1,3))
for(i in 1:length(transitions)){
  plot(etm.idm, tr.choice = transitions[i], conf.int = TRUE,
       lwd = 2, legend = FALSE, ylim = c(0, 1),
       xlim = c(0, 100), xlab = "Years",
       ci.fun = "cloglog")
}
par(mfrow=c(1,1))

summary(etm.idm)$"0 1"[, c("P", "lower", "upper")]

# landmark method for probability plots
time.points <- c(seq(3, 10, 1), 15)
landmark.etm <- lapply(time.points, function(start.time) {
    etm(my.icu.pneu, c("0", "1", "2"), tra.idm, "cens", start.time)
    })

# PAMMTOOLS ----

# data transformation
pamm.icu.pneu <- my.icu.pneu %>%
  mutate(status = ifelse(grepl("cens", transition), 0, 1)) %>%
  rename(tstart = entry
         , tstop = exit) %>%
  select(id
         , tstart
         , tstop
         , from
         , to
         , transition
         , status
         , age
         , sex)

my.pamm.icu.pneu_test <- pamm.icu.pneu %>% add_counterfactual_transitions()

cal_icu.pneu <- as_ped_multistate(
  data       = my.pamm.icu.pneu_test,
  formula    = Surv(tstart, tstop, status)~ age + sex + transition,
  transition = "transition",
  id         = "id",
  r = 0,
  timescale  = "calendar")

dim(pamm.icu.pneu)
dim(my.pamm.icu.pneu_test)
dim(cal_icu.pneu)

# predict hazards
# as factor: transition in cal_icu.pneu
# ggf in as_ped gleich als faktorvariable ausgeben.
pam <- mgcv::gam(
  formula = ped_status ~ s(tend, by=as.factor(transition)) + as.factor(transition) + s(age) + sex,
  data = cal_icu.pneu,
  family=poisson(),
  offset=offset)
summary(pam)

plot(pam, xlim = c(0,100), ylim = c(-20, 20), page=1)
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
