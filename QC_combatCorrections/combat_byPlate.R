## performs batch correction by plate (total of 7 plates)
## m values qc'd/preprocessed as detailed in "minfi_qc.R"
## note: this code is not designed to run directly 
### (files, paths, etc would need to be updated)

library(sva)
par(cex = 0.5)

mod_mat = model.matrix(~1, data=pData(gset.norm))
batch_pheno = gset.norm$Sample_Plate

M.batch <- ComBat(M.norm, batch=batch_pheno, mod_mat, par.prior = TRUE, prior.plots = FALSE)
## Found7batches

## Adjusting for0covariate(s) or covariate level(s)

## Standardizing Data across genes

## Fitting L/S model and finding priors

## Finding parametric adjustments

## Adjusting the Data

colnames(M.batch) <- sampleNames(gset.norm)

message("current dim of M.batch ", dim(M.batch)[1], " and ", dim(M.batch)[2], " samples, M normalized gset.")

## current dim of M.batch 439532 and 588 samples, M normalized gset.

