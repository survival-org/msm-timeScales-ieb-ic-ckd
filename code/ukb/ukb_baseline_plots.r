
# load packages ----
source("/wis37138/msm-timeScales-ieb-ic-ckd/code/helpers_ukb.r")

library(etm)
library(mvna)
library(kmi)
library(mgcv)
library(scales)
library(nnet)

# set directories ----
setwd("wis37138") # only necessary to enable plotting because I have no write permissions in "/"

dir_data <- "/wis37138/msm-timeScales-ieb-ic-ckd/data/ukb"
dir_figures <- "/wis37138/msm-timeScales-ieb-ic-ckd/results/ukb/msm/figures/export"


# 0>1 and 0>4 plot ----
plots_pam_stss_01 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_01.rds")))
plots_pam_stss_03 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_04.rds")))
plots_pam_mts_01 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_01.rds")))
plots_pam_mts_03 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_04.rds")))

# 1>2 and 1>4 plot ----
plots_pam_stss_age_12 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_age_12.rds")))
plots_pam_stss_age_13 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_age_14.rds")))
plots_pam_mts_age_12 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_age_12.rds")))
plots_pam_mts_age_13 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_age_14.rds")))

# 2>3 plot ----
plots_pam_stss_age_23 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_age_23.rds")))
plots_pam_mts_age_23 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_age_23.rds")))

# 2>4 plot ----
plots_pam_stss_age_24 <- readRDS(file.path(dir_data, paste0("pam_stss_plots_age_24.rds")))
plots_pam_mts_age_24 <- readRDS(file.path(dir_data, paste0("pam_mts_plots_age_24.rds")))