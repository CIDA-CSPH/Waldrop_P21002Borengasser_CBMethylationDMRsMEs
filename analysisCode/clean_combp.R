## we read in the results from combp to clean/format
## add some additional info

## Read in comb-p annotation results
outcome <- "adiposity_change_var"
# merged_for_bed, bed_combp ## see combp_setup.R
load(file = paste0("dataProcessed/outcome/outcome_combp_results.Rdata"))

bed_combp_results <- list()
for(i in 1:length(bed_combp)){
	bed_combp_results[[names(bed_combp)[i]]] <- fread(file = paste0("bedfiles/outcome_", names(bed_combp[i]), "BED_COMBP", ".bed_sorted.bed_forCombp.bed.combp_results.anno.hg19.bed"), 
																										header = TRUE,                            
																										sep = "\t")
}


## sort results by sidak p val
cb_cp_sort <- lapply(bed_combp_results, FUN = function(x){
	a <- x %>%
		arrange(z_sidak_p)
})
names(cb_cp_sort) <- names(bed_combp)

##FUNCTION methyl_dir_fun
## calculates direction of methylation
## number and % of positive (t_value > 0), negative t_value < 0
source(methyl_dir_fun.R)

cb_results <- list()
for(i in names(cb_cp_sort))
{
	# set dmp and dmr data to data.table
	dmp <- setDT(merged_for_bed[[i]])
	dmr <- setDT(cb_cp_sort[[i]])
	# set chrom as key for dmp.dt, faster binary search
	setkey(dmp, chrom, pos)
	setkey(dmr, "#chrom", start, end)
	
	myCluster <- makeCluster(15, # number of cores to use 
													 type = "FORK") # type of cluster
	doParallel::registerDoParallel(cl = myCluster)
	cb_results[[i]] <- foreach(j=1:nrow(dmr), .combine="rbind") %dopar%
		{
			methyl_dir_fun(x=j, dmr.dt=dmr, dmp.dt= dmp)
		}
	stopCluster(myCluster)
}

# save(cb_results, file = paste0("dataProcessed/outcome/COMBP_SORTED_RESULTS_outcome.RData"))
