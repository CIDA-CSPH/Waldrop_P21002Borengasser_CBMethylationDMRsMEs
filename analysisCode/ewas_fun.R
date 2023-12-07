# EWAS function

## Build function for model fitting in methylation, as predictor
# builds a linear regression model for each methylation probe,
# adjusting for covariates (covars), and selects the coefficient estimates, st error, t-value, and p-value to create a data frame of results for each linear regression.

# then, adjusts for multiple testing for all p-values from result data frame using Benjamini-Hochberg method

# cleans the results
# 1. coerces into data frame with probe, result data, and adjusted p-values
# 2. arrange data frame from smallest to largest adjusted p-value
# 3. convert data into numeric (first character, then numeric; p-value results can be non-numeric)

# arguments
# 1.beta: maternal clinical variable of interest, string
# 2. covars: input covariate values as argument
#   format covariate values as character string
#   specify any factor covariates in argument (ie use "as.factor(infant_sex)" in argument)
# 3.dat: clinical data (with cb data)
# 4.cpg: colnames of methylation data/probes
# 5.cluster: select the object to export to cluster (clin_dat_msc or clin_dat_plc), string

'%!in%' <- function(x,y)!('%in%'(x,y))
library(tidyverse)
library(data.table)
library(parallel)
library(devtools)
library(pbmcapply)

ewas_outcomes <- function(outcome, covars, dat, cpg, cluster){
	dat <- dat
	# col_vars <- c(cpg, outcome)
	# outcome <- dat[[outcome]]
	outcome <- dat[[outcome]]
	covars <- covars
	cpg <- cpg
	# define a local function, f, to perform a linear regression for predictor and extract the coefficient data from each regression into a resulting data frame. function will be called in a pbmclapply to run through all probes in cpg (msc)
	f <- local(function(x, outcome, dat){
		dat <- dat
		covars <- covars
		# linear regression of maternal clinical predictor, adjust for covars 
		if (length(covars) > 0)
		{
			# form <- as.formula(paste(outcome, " ~ ", x, " + ", paste(covars, collapse = "+")))
			form <- as.formula(paste("outcome ~ ", x, " + ", paste(covars, collapse = "+")))
		}
		else
		{
			form <- as.formula(paste("outcome ~ ", x))
		}
		
		fit1 <- lm(formula = form, data=dat, na.action = na.omit)
		fit.sum <- summary(fit1)
		cpg.dat <- fit.sum$coefficients[2, ] 
		return(cpg.dat)
	})
	
	# ## Set cluster
	cl <- parallel::makeCluster(20)
	parallel::clusterExport(cl, cluster)
	
	# form a data frame with linear regression coefficients for each probe (rows)
	rslts <- do.call("rbind", pbmclapply(cpg, FUN = function(x) f(x, outcome, dat), mc.cores =getOption("mc.cores", 20L)))
	
	## calculate FDR
	# Benjamini-Hochberg method
	p_adjust <- p.adjust(rslts[, 4], method = "BH")
	
	# create data frame with probe information, coefficient results from linear regressions, and adjusted p-values
	rslts <- as.data.frame(cbind(cpg, rslts, p_adjust))
	colnames(rslts)[1] <- "CpG"
	# arrange result data frame with adjusted p-values in ascending order
	rslts <- rslts %>% arrange(p_adjust)
	stopCluster(cl)
	return(rslts)
}