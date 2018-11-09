# PGC PTSD Epigenetics Working Group Quality Control and Analysis Pipelines

The following scripts can be used to run the quality control pipeline as detailed in [Ratanatharathorn et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.1002/ajmg.b.32568). 

Briefly, background Corrected beta-values, methylation signals, and detection p-values are loaded from iDATs in GenomeStudio and extracted into a txt file for QC to be performed in R using the following scripts:

1. 01_bkgdcor_QC_prep.R - extracts the beta-values, methylation signals, and detection p-values from the GenomeStudio output

2. 02_bkgdcor_QC_CpGassoc.R - Samples with probe detection call rates <90% and those with an average intensity value of either <50% of the experiment-wide sample mean or <2,000 arbitrary units (AU) are excluded. Data points with probe detection p-values >0.001 are set to missing, and CpG sites with missing data for >10% of samples are excluded from analysis.  

3. 03_bkgdcor_beta_BMIQ_bySample.R - Probes that cross hybridize between autosomes and sex chromosomes (Chen 2013) are removed and [Beta Mixture Quantile Normalization (BMIQ)](https://academic.oup.com/bioinformatics/article/29/2/189/204142) is run.

4. 04_bkgd_beta_ComBat_bySample.R - ComBat run to account for sources of technical variations.  Prior to ComBat correction, missing data is imputed using the nearest-neighbor method with default parameters.  Before imputation, a flag matrix is created with a true value in a particular cell indicating the presence of missing information for the particular cell. Beta values are then transformed into m-values. ComBat is run once to adjust for chip designation and a second time to adjust for the 12-point position designation. If the chips are not balanced, case designation and/or other relevant covariates (i.e. gender) are included. Following ComBat, m-values are converted back to beta values and missing values added back in using the flag matrix.

Two other scripts are used:

1. bkgdcor_beta_QC_popStrat_1bp.R - runs the [Barfield et al. (2014)](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.21789) code for estimating ancestry PCs using CpG sites within 1 bp of a SNP. PCs 2-4 are included in as covariates in the cross-sectional analysis.

2. FlowSortedBlood450K.R - runs the [minfi estimateCellCounts()](https://bioconductor.org/packages/devel/bioc/manuals/minfi/man/minfi.pdf) function to estimate cell types using the raw iDAT files.
