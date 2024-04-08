# details qc steps performed for cb methylation data
## Note this will not run directly; file paths must be updated

# Steps:
# Removal of samples with low median intensity (1 sample removed)
# Removal of probes with a detection p-value greater than 0.01 in more than 10% of samples (699 probes removed)
# Removal of samples with a detection p-value greater than 0.01 in more than 1% of probes (0 samples removed)
# Removal of probes with a bead count <3 in at least 5% of samples (660 probes removed)
# Non-CG probes were assessed but were not removed (3,091 non-CG probes kept)
# Non-mapping probes were assessed (0 non-mapping probes)
# Removal of probes with SNPs at the CpG interrogation and/or at the single nucleotide extension for any minor allele frequency (17,272 probes removed)
# Removal of cross-reactive probes (27,349 probes removed. The list of known cross-reactive probes was identified in the paper by Chen et al. (https://www.tandfonline.com/doi/full/10.4161/epi.23470)
# Removal of samples with mismatched sex between the clinical data and predicted sex through the minfi getSex function (6 samples removed)
# Allosomes and autosomes were noted. The CpG names for both autosomes (429,246 probes) and allosomes (10,286 probes) can be found on the server for reference (“autosomes.txt” and “allosomes.txt”).

rm(list = ls())
library(knitr)
options(stringsAsFactors = F)
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(eval = TRUE)
knitr::opts_chunk$set(message = TRUE)
knitr::opts_chunk$set(warning = FALSE)
knitr::opts_chunk$set(fig.show = 'hold')
knitr::opts_chunk$set(fig.align = 'center')
knitr::opts_chunk$set(out.width = '100%')
knitr::opts_chunk$set(out.height = '100%')
library(rmarkdown)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(FlowSorted.Blood.450k)
library(wateRmelon)
library(sva)
library(plyr)
library(data.table)
library(outliers)
library(RColorBrewer)

'%!in%' <- function(x,y)!('%in%'(x,y))

# read in data from each plate using  read.metharray.sheet()
##Load data
# plate 1
base1dir <- paste0(dir,"Plate_1_04232015/")
targets1 <- read.metharray.sheet(base1dir)

# plate 2
base2dir <- paste0(dir,"plate_2_04222015/")
targets2 <- read.metharray.sheet(base2dir)

# plate 3
base3dir <- paste0(dir,"plate_3_04292015/")
targets3 <- read.metharray.sheet(base3dir)

# plate 4
base4dir <- paste0(dir,"plate4_04302015/")
targets4 <- read.metharray.sheet(base4dir)

# plate 5
base5dir <- paste0(dir,"plate5_05052015/")
targets5 <- read.metharray.sheet(base5dir)

# plate 6
base6dir <- paste0(dir,"plate6_05062015/")
targets6 <- read.metharray.sheet(base6dir)

# plate 7
base7dir <- paste0(dir,"plate7_05082015/")
targets7 <- read.metharray.sheet(base7dir)

targets <- rbind(targets1, targets2, targets3, targets4, targets5, targets6, targets7) # n=600

rgSet <- read.metharray.exp(targets=targets, extended=T)

# update pid to reflect offspring
rgSet[[1]] <- rgSet[[1]] + 10000
# replace mis-labeled 50568 pid to 40568
rgSet[[1]][rgSet[[1]]==50568] <- 40568
sampleNames(rgSet) = rgSet[[1]]

rgSet

## class: RGChannelSetExtended 
## dim: 622399 600 
## metadata(0):
## assays(5): Green Red GreenSD RedSD NBeads
## rownames(622399): 10600313 10600322 ... 74810490 74810492
## rowData names(0):
## colnames(600): 20850 20686 ... 21322 21396
## colData names(9): Sample_Name Sample_Well ... Basename filenames
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19

message("initial rgSet: ", dim(rgSet)[2], " samples")

## initial rgSet: 600 samples

# save(rgSet, file=paste0(dir_ext, "rgSet.Rdata"))

# qc report, includes qc plots from minfi
qcReport(rgSet, sampNames=targets$Sample_Name, sampGroups=targets$Sample.Group, 
				 pdf=paste0(dir_ext, "qcReport.pdf"))

## Warning: NON-POLYMORPHIC probes outside plot range

getManifest(rgSet)

## IlluminaMethylationManifest object
## Annotation
##   array: IlluminaHumanMethylation450k
## Number of type I probes: 135476 
## Number of type II probes: 350036 
## Number of control probes: 850 
## Number of SNP type I probes: 25 
## Number of SNP type II probes: 40

clin_dat <- fread("/home/smiharry/LEAD/dataraw_HS/k_boyle/clinical/boyle190419.csv") # 651
clin_dat$PID <- clin_dat$PID + 10000
## check pid in ECHOIDs file (provided )
revoke <- read.csv("/home/niemiecs/ECHO_Aim2/Methylation/cb/cb.revoked.pid.csv", header=T) # 5 pid that revoked consent
# variables
#names(clin_dat)
nrow(clin_dat)

## [1] 651

message("Clin_dat Dim: ", dim(clin_dat)[[1]], " samples, ", dim(clin_dat)[[2]], " variables")

## Clin_dat Dim: 651 samples, 151 variables

clindat2 <- clin_dat[match(rgSet$Sample_Name, clin_dat$PID),]
# rm revoked
clindat2 <- clindat2[which(clindat2$PID %!in% revoke$revokedid),] # n = 595

rgSet <- rgSet[, rgSet$Sample_Name %in% clindat2$PID] # n = 595
rgSet

## class: RGChannelSetExtended 
## dim: 622399 595 
## metadata(0):
## assays(5): Green Red GreenSD RedSD NBeads
## rownames(622399): 10600313 10600322 ... 74810490 74810492
## rowData names(0):
## colnames(595): 20850 20686 ... 21322 21396
## colData names(9): Sample_Name Sample_Well ... Basename filenames
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19

stopifnot(all(clindat2$PID == rgSet$Sample_Name))


head(pData(rgSet))

## DataFrame with 6 rows and 9 columns
##       Sample_Name Sample_Well Sample_Plate Sample.Group     Pool_ID       Array
##         <numeric> <character>  <character>    <logical> <character> <character>
## 20850       20850         A01      Plate 1           NA          NA      R01C01
## 20686       20686         B01      Plate 1           NA          NA      R02C01
## 20667       20667         C01      Plate 1           NA          NA      R03C01
## 20676       20676         D01      Plate 1           NA          NA      R04C01
## 20688       20688         E01      Plate 1           NA          NA      R05C01
## 20666       20666         F01      Plate 1           NA          NA      R06C01
##             Slide
##       <character>
## 20850  3999120006
## 20686  3999120006
## 20667  3999120006
## 20676  3999120006
## 20688  3999120006
## 20666  3999120006
##                                                                              Basename
##   (printed details on file paths rm here)                                                                        <character>

#### clindat_pheno = clindat2 %>% select(c("PID", "infant_sex"))
#### pData(rgSet) = merge(pData(rgSet), clindat_pheno, by.x="Sample_Name", by.y="PID")
# identical(pData(rgSet)$Sample_Name, clindat2$PID)

## [1] TRUE

pData(rgSet)$infant_sex <- ifelse(clindat2$infant_sex == 1, "F", "M")
message("Clin_dat2 Dim: ", dim(clindat2)[[1]], " samples, ", dim(clindat2)[[2]], " variables")

## Clin_dat2 Dim: 595 samples, 151 variables
message("rgSet: ", dim(rgSet)[2], " samples")

## rgSet: 595 samples

cell.counts <- estimateCellCounts(rgSet=rgSet, meanPlot=T, compositeCellType="CordBlood", cellTypes=c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC"))

## Loading required package: FlowSorted.CordBlood.450k

## [estimateCellCounts] Combining user data with reference (flow sorted) data.

## [estimateCellCounts] Processing user and reference data together.

## [estimateCellCounts] Picking probes for composition estimation.

## [estimateCellCounts] Estimating composition.

cell.counts <- data.table::as.data.table(cell.counts, keep.rownames = "pid")

cell.counts[1:10,]

##       pid      CD8T      CD4T          NK      Bcell       Mono        Gran
##  1: 20850 0.2857422 0.2081771 0.036354756 0.13081950 0.11253353 0.004560078
##  2: 20686 0.2037473 0.2159390 0.067477957 0.09123277 0.13571422 0.261161315
##  3: 20667 0.1874163 0.1602650 0.008294263 0.10761988 0.05833561 0.429726421
##  4: 20676 0.1287379 0.1104240 0.000000000 0.05479837 0.12672549 0.537805407
##  5: 20688 0.2415234 0.2206213 0.124569533 0.08317922 0.08971951 0.196096897
##  6: 20666 0.1787222 0.1207732 0.037341128 0.09104510 0.14339290 0.362973184
##  7: 20672 0.2066902 0.2440690 0.169225774 0.19128022 0.11389126 0.070322111
##  8: 20674 0.3292179 0.3149130 0.000000000 0.12726769 0.10486025 0.104103314
##  9: 20892 0.1069909 0.1508287 0.000000000 0.08655135 0.06578939 0.548665221
## 10: 20654 0.2498828 0.3278061 0.053893449 0.19503222 0.10617849 0.048829639
##           nRBC
##  1: 0.24336346
##  2: 0.06558382
##  3: 0.10140953
##  4: 0.06655826
##  5: 0.05771854
##  6: 0.10552262
##  7: 0.04840490
##  8: 0.05915894
##  9: 0.09526054
## 10: 0.06949887

# fwrite(cell.counts, file=paste0(dir_ext, "/dataprocessed/cell.counts.csv"))

mset <- preprocessRaw(rgSet = rgSet)
mset

## class: MethylSet 
## dim: 485512 595 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(485512): cg00050873 cg00212031 ... ch.22.47579720R
##   ch.22.48274842R
## rowData names(0):
## colnames(595): 20850 20686 ... 21322 21396
## colData names(10): Sample_Name Sample_Well ... filenames infant_sex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

qc <- getQC(mset)

# outliers: [Q1 - 3*IQR, Q3 + 3*IQR]
# assess for far outliers for methylated median intensity
# upper outliers, 0
#message("Methylated, Far-Out Upper Outliers: ", rownames(qc)[qc$mMed > (quantile(qc$mMed, 0.75) + 3*IQR(qc$mMed))])
# lower outliers, n=1 "20772"
# rownames(qc)[qc$mMed < (quantile(qc$mMed, 0.25) - 3*IQR(qc$mMed))]
badid <- rownames(qc)[qc$mMed < (quantile(qc$mMed, 0.25) - 3*IQR(qc$mMed))]
message("Methylated, Far-Out Lower Outliers: ", rownames(qc)[qc$mMed < (quantile(qc$mMed, 0.25) - 3*IQR(qc$mMed))])

## Methylated, Far-Out Lower Outliers: 20772

# assess for far outliers for unmethylated median intensity
# message("Unmethylated, Far-Out Upper Outliers: ", rownames(qc)[qc$mMed > (quantile(qc$mMed, 0.75) + 3*IQR(qc$mMed))])
# message("Unmethylated, Far-Out Lower Outliers: ", rownames(qc)[qc$uMed < (quantile(qc$uMed, 0.25) - 3*IQR(qc$uMed))])

# confirm with grubbs, normality assumption satisfied
library(outliers)
#methylated
grubbs.test(qc$mMed) # significant, lowest value is outlier

## 
##  Grubbs test for one outlier
## 
## data:  qc$mMed
## G = 5.02757, U = 0.95738, p-value = 0.0001121
## alternative hypothesis: lowest value 11.7427304765192 is an outlier


#unmethylated
grubbs.test(qc$uMed) # insignificant, no outliers

## 
##  Grubbs test for one outlier
## 
## data:  qc$uMed
## G = 3.20474, U = 0.98268, p-value = 0.3849
## alternative hypothesis: lowest value 11.5584207132687 is an outlier

message("Bad sample: ", badid)

## Bad sample: 20772

write.table(badid, file=paste0(dir_ext, "badid.txt"), quote=F, row.names=F, col.names=F)
message("Number of bad samples labeled: ", length(badid))

## Number of bad samples labeled: 1

mset <- mset[,which(mset$Sample_Name %!in% badid)]

message("Current dim of mset: ", dim(mset)[1], " probes and ", dim(mset)[2], " samples, following removal of probes with low beadcount")

## Current dim of mset: 485512 probes and 594 samples, following removal of probes with low beadcount

## Investigate bad probes
detP <- detectionP(rgSet = rgSet)
detPcut <- 0.01

failed <- detP > detPcut
numfail <- colMeans(failed)
x <- seq(1, length(numfail), 1)

x[numfail==max(numfail)]

## [1] 232

removeDetP <- 0.1
badProbes <- rowMeans(failed) > removeDetP
badProbeNames1 <- names(badProbes[which(badProbes==T)])

mset.f <- mset[!badProbes,]
mset.f

## class: MethylSet 
## dim: 484813 594 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(484813): cg00455876 cg01707559 ... ch.22.47579720R
##   ch.22.48274842R
## rowData names(0):
## colnames(594): 20850 20686 ... 21322 21396
## colData names(10): Sample_Name Sample_Well ... filenames infant_sex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

message("Probes with detection p-value > 0.01 in more than 10% of samples: ", dim(mset)[1] - dim(mset.f)[1], " probes")

## Probes with detection p-value > 0.01 in more than 10% of samples: 699 probes

mset = mset.f

message("Current dim of mset, probes: ", dim(mset)[1], " and ", dim(mset)[2], " samples, following removal of probes with low p-value")

## Current dim of mset, probes: 484813 and 594 samples, following removal of probes with low p-value


samps_numfail <- rowMeans(failed)

x[samps_numfail==max(samps_numfail)]

## [1] 205120

samp_removeDetP <- 0.01
badSamps <- colMeans(failed) > removeDetP
message("Number of bad samples: ", sum(badSamps==T))

## Number of bad samples: 0

sum(badSamps==T)

## [1] 0

message("Samples with detection p-value > 0.01 in more than 1% of samples: ", sum(badSamps==T))

## Samples with detection p-value > 0.01 in more than 1% of samples: 0

message("Current dim of mset: ", dim(mset)[1], " probes and ", dim(mset)[2], " samples, following remove samples with low p-value")

## Current dim of mset: 484813 probes and 594 samples, following remove samples with low p-value


beadCutoff = 0.05
bc <- beadcount(rgSet)
quantile(bc, na.rm=T)

##   0%  25%  50%  75% 100% 
##    3   11   13   16  255

bc2 = bc[rowSums(is.na(bc)) < beadCutoff * (ncol(bc)), ]
badProbeNames2 <- row.names(bc[rowSums(is.na(bc)) >= beadCutoff*(ncol(bc)), ])
mset.f2 = mset[featureNames(mset) %in% row.names(bc2), ]
mset.f2

## class: MethylSet 
## dim: 484153 594 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(484153): cg00455876 cg01707559 ... ch.22.47579720R
##   ch.22.48274842R
## rowData names(0):
## colnames(594): 20850 20686 ... 21322 21396
## colData names(10): Sample_Name Sample_Well ... filenames infant_sex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

message("Filtering probes with a beadcount <3 in at least ", beadCutoff*100, "% of samples, removed ", dim(mset)[1] - dim(mset.f2)[1], " from the analysis.")

## Filtering probes with a beadcount <3 in at least 5% of samples, removed 660 from the analysis.

mset = mset.f2

message("Current dim of mset: ", dim(mset)[1], " probes and ", dim(mset)[2], " samples, following removal of probes with low beadcount")

## Current dim of mset: 484153 probes and 594 samples, following removal of probes with low beadcount

mset.cg <- dropMethylationLoci(mset, dropCH=T)
mset.cg

## class: MethylSet 
## dim: 481065 594 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(481065): cg00455876 cg01707559 ... cg27662611 cg27665648
## rowData names(0):
## colnames(594): 20850 20686 ... 21322 21396
## colData names(10): Sample_Name Sample_Well ... filenames infant_sex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

message("Non-CG Probes: ", dim(mset)[1] - dim(mset.cg)[1], ". These were kept.")

## Non-CG Probes: 3088. These were kept.


gset <- mapToGenome(mset)
gset

## class: GenomicMethylSet 
## dim: 484153 594 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(484153): cg13869341 cg14008030 ... cg01757887 cg21106100
## rowData names(0):
## colnames(594): 20850 20686 ... 21322 21396
## colData names(10): Sample_Name Sample_Well ... filenames infant_sex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

message("Dim of gset: ", dim(gset)[1], " probes and ", dim(gset)[2], " samples")

## Dim of gset: 484153 probes and 594 samples

annotation <- getAnnotation(gset, dropNonMapping=F)
names(annotation)

##  [1] "chr"                      "pos"                     
##  [3] "strand"                   "Name"                    
##  [5] "AddressA"                 "AddressB"                
##  [7] "ProbeSeqA"                "ProbeSeqB"               
##  [9] "Type"                     "NextBase"                
## [11] "Color"                    "Probe_rs"                
## [13] "Probe_maf"                "CpG_rs"                  
## [15] "CpG_maf"                  "SBE_rs"                  
## [17] "SBE_maf"                  "Islands_Name"            
## [19] "Relation_to_Island"       "Forward_Sequence"        
## [21] "SourceSeq"                "Random_Loci"             
## [23] "Methyl27_Loci"            "UCSC_RefGene_Name"       
## [25] "UCSC_RefGene_Accession"   "UCSC_RefGene_Group"      
## [27] "Phantom"                  "DMR"                     
## [29] "Enhancer"                 "HMM_Island"              
## [31] "Regulatory_Feature_Name"  "Regulatory_Feature_Group"
## [33] "DHS"

table(annotation$chr)

## 
##  chr1 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19  chr2 chr20 
## 46771 24341 28732 24483 12264 15039 15235 21926 27827  5911 25456 34739 10349 
## chr21 chr22  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chrX  chrY 
##  4237  8522 25106 20430 24276 36484 29954 20913  9849 11206   103

dim(annotation)

## [1] 484153     33

annotation2 <- getAnnotation(gset, dropNonMapping=T)
dim(annotation2)

## [1] 484153     33

message("non-mapping Probes: ", dim(annotation)[1] - dim(annotation2)[1], ". These were kept.")

## non-mapping Probes: 0. These were kept.#0

message("Current dim of gset: ", dim(gset)[1], " probes and ", dim(gset)[2], " samples, following non-mapping")

## Current dim of gset: 484153 probes and 594 samples, following non-mapping


snps <- getSnpInfo(gset)
head(snps, 10)

## DataFrame with 10 rows and 6 columns
##               Probe_rs Probe_maf      CpG_rs   CpG_maf      SBE_rs   SBE_maf
##            <character> <numeric> <character> <numeric> <character> <numeric>
## cg13869341          NA        NA          NA        NA          NA        NA
## cg14008030          NA        NA          NA        NA          NA        NA
## cg12045430          NA        NA          NA        NA          NA        NA
## cg20826792          NA        NA          NA        NA          NA        NA
## cg00381604          NA        NA          NA        NA          NA        NA
## cg20253340          NA        NA          NA        NA          NA        NA
## cg21870274          NA        NA          NA        NA          NA        NA
## cg03130891  rs77418980  0.305556          NA        NA          NA        NA
## cg24335620 rs147502335  0.012800          NA        NA          NA        NA
## cg16162899          NA        NA          NA        NA          NA        NA

dim(snps)

## [1] 484153      6

gset <- addSnpInfo(gset)
head(granges(gset))

## GRanges object with 6 ranges and 6 metadata columns:
##              seqnames    ranges strand |    Probe_rs Probe_maf      CpG_rs
##                 <Rle> <IRanges>  <Rle> | <character> <numeric> <character>
##   cg13869341     chr1     15865      * |        <NA>        NA        <NA>
##   cg14008030     chr1     18827      * |        <NA>        NA        <NA>
##   cg12045430     chr1     29407      * |        <NA>        NA        <NA>
##   cg20826792     chr1     29425      * |        <NA>        NA        <NA>
##   cg00381604     chr1     29435      * |        <NA>        NA        <NA>
##   cg20253340     chr1     68849      * |        <NA>        NA        <NA>
##                CpG_maf      SBE_rs   SBE_maf
##              <numeric> <character> <numeric>
##   cg13869341        NA        <NA>        NA
##   cg14008030        NA        <NA>        NA
##   cg12045430        NA        <NA>        NA
##   cg20826792        NA        <NA>        NA
##   cg00381604        NA        <NA>        NA
##   cg20253340        NA        <NA>        NA
##   -------
##   seqinfo: 24 sequences from hg19 genome; no seqlengths

getAnnotationObject(gset)

## IlluminaMethylationAnnotation object
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12
##   genomeBuild: hg19
## Available annotation
##   Islands.UCSC
##   Locations
##   Manifest
##   Other
##   SNPs.132CommonSingle
##   SNPs.135CommonSingle
##   SNPs.137CommonSingle
##   SNPs.138CommonSingle
##   SNPs.141CommonSingle
##   SNPs.142CommonSingle
##   SNPs.144CommonSingle
##   SNPs.146CommonSingle
##   SNPs.147CommonSingle
##   SNPs.Illumina
## Defaults
##   Locations
##   Manifest
##   SNPs.137CommonSingle
##   Islands.UCSC
##   Other

gset.f2 <- dropLociWithSnps(gset, snps=c("SBE", "CpG"), maf=0)
badProbeNames3 <- featureNames(gset)[!featureNames(gset) %in% featureNames(gset.f2)]

message("There were ", dim(gset)[1] - dim(gset.f2)[1], " probes with snps.")

## There were 17272 probes with snps.

gset <- gset.f2
dim(gset)

## [1] 466881    594

message("Current dim of gset: ", dim(gset)[1], " probes and ", dim(gset)[2], " samples, following removal of probes with snps")

## Current dim of gset: 466881 probes and 594 samples, following removal of probes with snps

# cross reactive probe file
xReactiveProbes <- read.csv("/home/datasets/cross_reactive/non_specific_cg.csv")
dim(xReactiveProbes)[1]

## [1] 29233

keep <- !(featureNames(gset) %in% xReactiveProbes$TargetID)
table(keep)

## keep
##  FALSE   TRUE 
##  27349 439532

xReactive_data <- gset[which(featureNames(gset) %in% xReactiveProbes$TargetID),]
message("There were ", dim(xReactive_data)[1], " cross-reactive probes.")

## There were 27349 cross-reactive probes.

gset.f3 <- gset[keep, ]
dim(gset.f3)

## [1] 439532    594

gset <- gset.f3

message("Current dim of gset: ", dim(gset)[1], " probes and ", dim(gset)[2], " samples, following removal of cross-reactive probes")

## Current dim of gset: 439532 probes and 594 samples, following removal of cross-reactive probes

# paper by Chen et al. (<https://www.tandfonline.com/doi/full/10.4161/epi.23470>). 


beta.raw <- getBeta(gset)
M.raw <- getM(gset)


snames <- sampleNames(gset)

## Get median intensities for X and Y chroms
gset.sex <- getSex(gset) 
## Add these columns to the phenotype data for the gset
pData(gset) <- cbind(pData(gset), gset.sex)
phendata <- pData(gset)
## create a separate phenotype dataset

##Compare predicted sex with clinical sex data
badsex <- phendata[which(phendata$predictedSex != phendata$infant_sex), 1] # n = 6
write.table(badsex, file=paste0(dir_ext, "badsex.txt"), quote = F, row.names = F, col.names = F)
## remove these samples from dataset 

pheno.clean <- phendata[which(phendata$Sample_Name %!in% badsex), ]
gset.clean <- gset[, gset$Sample_Name %!in% badsex] # n = 589

message("There were ", dim(gset)[2] - dim(gset.clean)[2], " samples with mismatched sex.")

## There were 6 samples with mismatched sex.

gset.clean

## class: GenomicMethylSet 
## dim: 439532 588 
## metadata(0):
## assays(2): Meth Unmeth
## rownames(439532): cg13869341 cg24669183 ... cg25918849 cg21106100
## rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
## colnames(588): 20850 20686 ... 21322 21396
## colData names(13): Sample_Name Sample_Well ... yMed predictedSex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

gset = gset.clean
message("Current dim of gset ", dim(gset)[1], " and ", dim(gset)[2], " samples, removal of sex mismatch.")

## Current dim of gset 439532 and 588 samples, removal of sex mismatch.

message("PIDs with mismatched sex ", badsex[1], ", ", badsex[2], ", ",badsex[3], ", ",badsex[4], ", ",badsex[5], ", ",badsex[6])

## PIDs with mismatched sex 20661, 20696, 20723, 20791, 21408, 40568


# probes not on sex chromosomes
autosomes <- annotation[!annotation$chr %in% c("chrX", "chrY"),]
# probes on sex chromosomes 
allosomes <- annotation[annotation$chr %in% c("chrX", "chrY"),]

# filter out allosomes
gset.auto <- gset[featureNames(gset) %in% row.names(autosomes), ]

message("Autosomes: ", dim(gset.auto)[1], " within the dataset.")

## Autosomes: 429246 within the dataset.

message("Allosomes: ", dim(gset)[1] - dim(gset.auto)[1], "  within the dataset.")

## Allosomes: 10286  within the dataset.

## 10,286
message("Dim of autosomes gset ", dim(gset.auto)[1], " and ", dim(gset.auto)[2], " samples, following removal of allosomes.")

## Dim of autosomes gset 429246 and 588 samples, following removal of allosomes.

# write files
autosomes_data = data.frame(rownames(gset.auto))
colnames(autosomes_data) <- "cpg"
# write.table(autosomes_data, file=paste0(dir_ext, "autosomes.txt"), row.names = F, col.names =T)

allosomes_data = data.frame(rownames(gset)[which(rownames(allosomes) %in% rownames(gset))])
colnames(allosomes_data) <- "cpg"
# write.table(allosomes_data, file=paste0(dir_ext, "allosomes.txt"), row.names = F, col.names = T)

# raw beta
gset.norm <- preprocessQuantile(gset, removeBadSamples = T) # 439,532

## [preprocessQuantile] Mapping to genome.

## [preprocessQuantile] Fixing outliers.

## [preprocessQuantile] Quantile normalizing.

gset.norm

## class: GenomicRatioSet 
## dim: 439532 588 
## metadata(0):
## assays(2): M CN
## rownames(439532): cg13869341 cg24669183 ... cg25918849 cg21106100
## rowData names(6): Probe_rs Probe_maf ... SBE_rs SBE_maf
## colnames(588): 20850 20686 ... 21322 21396
## colData names(13): Sample_Name Sample_Well ... yMed predictedSex
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
## Preprocessing
##   Method: Raw (no normalization or bg correction)
##   minfi version: 1.34.0
##   Manifest version: 0.4.0

message("Current dim of gset.norm ", dim(gset.norm)[1], " and ", dim(gset.norm)[2], " samples, normalized gset.")

## Current dim of gset.norm 439532 and 588 samples, normalized gset.

## dim: 439,532, 589

# beta values
beta.norm <- getBeta(gset.norm)
# identical(colnames(beta.norm), sampleNames(gset))

message("Current dim of beta.norm ", dim(beta.norm)[1], " and ", dim(beta.norm)[2], " samples, beta normalized gset.")

## Current dim of beta.norm 439532 and 588 samples, beta normalized gset.

## dim: 439,532, 589

# M values
M.norm <- getM(gset.norm)
# identical(colnames(M.norm), sampleNames(gset))
# range(M.norm)
message("Current dim of M.norm ", dim(M.norm)[1], " and ", dim(M.norm)[2], " samples, M normalized gset.")

## Current dim of M.norm 439532 and 588 samples, M normalized gset.

## 439532 and 589 sample

# filter to only autosomes
gset.norm.auto <- gset.norm[featureNames(gset.norm) %in% row.names(autosomes),] # 429,246 (autosomes), difference of 10,286 (allosomes)
beta.norm.auto <- getBeta(gset.norm.auto)

M.norm.auto <- getM(gset.norm.auto)

message("Current dim of M.norm.auto ", dim(M.norm.auto)[1], " and ", dim(M.norm)[2], ", M normalized, autosomes only gset.")

## Current dim of M.norm.auto 429246 and 588, M normalized, autosomes only gset.

