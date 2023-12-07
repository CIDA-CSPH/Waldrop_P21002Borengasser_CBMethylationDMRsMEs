## this code details how ewas results were generated

## libraries
options(stringsAsFactors = FALSE)

'%!in%' <- function(x,y)!('%in%'(x,y))
library(tidyverse)
library(data.table)
library(parallel)
library(devtools)
library(pbmcapply)
library(labelled)

library(gridExtra)
library(RColorBrewer)
library(kableExtra)
library(cowplot)

# function to perform ewas
source("ewas_fun.R")

# cb : data table of methylation (M values, post-QC, post-batch correction)
# cb_clin: clinical data for subjects

# outcomes: change in adiposity by each time point (3 time points)
# outcome 1: adiposity IPV4 - adiposity IPV3 ("change_ipv4_ipv3")
# outcome 2: adiposity HSII - adiposity IPV3 ("change_hsII_ipv3")
# outcome 3: adiposity HSII - adiposity IPV4 ("change_hsII_ipv4")
## this was performed for each of these outcomes
# adiposity_change_var: captures change in adiposity measure used for model

# data
cb_cpgs <- colnames(cb)[which(colnames(cb) %!in% "pid")]

##################################################
# MODEL 1
##################################################
mod1_covars_vars <- c("race", "ga_at_birth_weeks", "infant_sex")
mod1_cb_fit <- ewas_outcomes(outcome = "adiposity_change_var", 
																			 covars = mod1_covars_vars, 
																			 dat = cb_clin,
																			 cpg = cb_cpgs, 
																			 cluster = "clin_dat")

##################################################
# MODEL 2
##################################################
mod2_cb_fit_ipv3 <- c("race", "ga_at_birth_weeks", "infant_sex", "pre_preg_bmi", "Gravidity", "gestsmoking")
mod2_cb_fit <- ewas_outcomes(outcome = "adiposity_change_var", 
																			 covars = mod2_cb_fit_ipv3, 
																			 dat = cb_clin,
																			 cpg = cb_cpgs, 
																			 cluster = "clin_dat")

##################################################
# MODEL 3
##################################################
mod3_covars_vars <- c("race", "ga_at_birth_weeks", "infant_sex", "pre_preg_bmi", "Gravidity", "gestsmoking", "DeliveryMode", "feedingType")
mod3_cb_fit <- ewas_outcomes(outcome = "adiposity_change_var", 
																			 covars = mod3_covars_vars, 
																			 dat = cb_clin,
																			 cpg = cb_cpgs, 
																			 cluster = "clin_dat")

change_cb_fit_list <- list(mod1_cb_fit, mod2_cb_fit, mod3_cb_fit)

names(change_cb_fit_list) <- c("mod1", "mod2", "mod3")
print(paste0("names of cb_fit_list_nogdm:", names(change_cb_fit_list)))
## Check that all tests were performed
sapply(change_cb_fit_list, function(x) length(cb_cpgs) == dim(x)[1]) ## Should all be TRUE

# save results, will use for combp_setup.R
save(change_cb_fit_list, file = paste0("ewas_results.Rdata"))

