library(batchtools)
library(checkmate)
library(foreach)
library(doParallel)

# setup ----

setwd("nvmetmp/wis37138/msm-timeScales-ieb-ic-ckd")
source("code/simulations/helpers_sim.r")
source("code/simulations/helpers_ic.r")
