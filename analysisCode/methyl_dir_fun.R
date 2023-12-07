##FUNCTION methyl_dir_fun
## calculates direction of methylation
## number and % of positive (t_value > 0), negative t_value < 0

# prior to call:
# set dmp and dmr data to data.table
## dmp <- setDT(dmp)
## dmr <- setDT(dmr)
# set chrom as key for dmp.dt, faster binary search
## setkey(dmp.dt, chrom)
## setkey(dmr.dt, "#chrom")

# function definition
# input
## x: row index of dmr data
## dmp: data.table
methyl_dir_fun <- function(x, dmr.dt, dmp.dt)
{
	# set dmr data
	chr <- dmr.dt$"#chrom"[[x]]
	start <-  dmr.dt$start[[x]]
	end <- dmr.dt$end[[x]]
	
	tmp <- dmp.dt[chr][pos %between% c(start, end)][, .(.N), keyby=.(pos=t_value > 0)]
	pos_n <- ifelse(length(tmp[pos==TRUE, N] > 0), tmp[pos==TRUE, N], 0)
	pos_proportion <- (pos_n/(sum(tmp$N)))
	neg_n <- ifelse(length(tmp[pos==FALSE, N] > 0), tmp[pos==FALSE, N], 0)
	neg_proportion <- (neg_n/(sum(tmp$N)))
	tmp.df <- cbind(dmr.dt[x,], pos_n, pos_proportion, neg_n, neg_proportion)
	return(tmp.df)
}
