# Waldrop_P21002Borengasser_CBMethylationDMRsMEs
This repository documents the code used to analyze the data for "Cord Blood DNA Methylation of Differentially Methylated Regions Associates with Infant and Child adiposity”

Methylation data underwent QC/preprocessing with the minfi pipeline, detailed in minfi_qc.R.

Combat was used for batch correction by sample plate (7 plates total), detailed in combat_byPlate.R

Methylation data was batch corrected using combat prior to analysis. We show the adequacy of this correction with MDS plots of probes by plate before and after correction in QC_combatCorrections/

Code used for analysis are detailed under analysisCode/

This code is generalized here, but was run for each of the outcomes of interest (changes in adiposity across the 3 timepoints, pairwise)

Analysis proceeded in the following order:
* ewas_code.R (includes function defined in ewas_fun.R)
* combp_setup.R
* combp_analyses.sh
* clean_combp.R (includes function defined in methyl_dir_fun.R)

Session information for these analyses are detailed in analysis_sessionInfo.md


