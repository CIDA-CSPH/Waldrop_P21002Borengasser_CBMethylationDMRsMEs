## this code details the compb analysis
## ewas results were loaded, see: ewas_code.R for code that details how these were calculated
## this analysis was performed a total of 3 times, once for each of the outcomes


# libraries for analysis
options(stringsAsFactors = FALSE)
'%!in%' <- function(x,y)!('%in%'(x,y))
library(tidyverse)
library(data.table)
library(parallel)
library(devtools)
library(pbmcapply)
library(gridExtra)
library(RColorBrewer)
library(kableExtra)
library(labelled)
# annotation for illumina 450k array
# cord blood run on 450k
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)


## Load hg19 annotation data
# illumina 450k
data("IlluminaHumanMethylation450kanno.ilmn12.hg19")
annotation.table <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
dim(annotation.table) #485512     33

## Bind CpG names with refSeq gene Names
# name refers to CpG
mapping <- data.frame("Name" = annotation.table$Name,
											"RefSeq" = annotation.table$UCSC_RefGene_Name,
											"chrom" = annotation.table$chr,
											"pos" = annotation.table$pos, 
											"strand" = annotation.table$strand)
dim(mapping) #485512      5


### DMR analysis with comb-p
## select outcome/model
## load ewas results: ewas_results, generated from ewas_code.R

## Merge mapping df with EWAS results, for analysis
merged_for_bed <- lapply(ewas_results, FUN = function(x){
	a <- merge(mapping, x, by.x = "Name", by.y = "CpG")
	colnames(a) <- c(names(a)[1:6], "std_error", "t_value", "rawp", "p_adjust")
	return(a)
})
names(merged_for_bed)
# "mod1" "mod2" "mod3"
names(merged_for_bed[[1]])
# [1] "Name"      "RefSeq"    "chrom"     "pos"       "strand"    "Estimate" 
# [7] "std_error" "t_value"   "rawp"      "p_adjust" 


## Generate bed files formatted for comb-p
# sorted by chromosome position
# total sample
bed_combp <- lapply(merged_for_bed, FUN= function(x){
	chrOrder <- c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
	a <- x %>%
		mutate(end = pos + 51,
					 chrom = factor(chrom, levels = chrOrder)) %>%
		select(chrom, pos, end, rawp) %>%
		mutate(pos=as.integer(pos), end=as.integer(end)) %>%
		arrange(chrom, pos)
})
names(bed_combp[[1]])
# [1] "chrom" "pos"   "end"   "rawp" 
dim(bed_combp[[1]])
# [1] 439532      4

## Write out BED files with no column name to run sortBed
for(i in 1:length(bed_combp)){
	fwrite(bed_combp[[i]], file = paste0("bedfiles/outcome_", names(bed_combp[i]), "BED_COMBP.bed"), 
				 row.names = FALSE,
				 col.names = FALSE,
				 sep = "\t")
}
