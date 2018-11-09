# PGC PTSD Epigenetics Working Group Quality Control and Analysis Pipelines

The following scripts can be used to run the quality control pipeline as detailed in [Ratanatharathorn et al. (2017)](https://onlinelibrary.wiley.com/doi/abs/10.1002/ajmg.b.32568). 

Briefly, background Corrected beta-values, methylation signals, and detection p-values are loaded from iDATs in GenomeStudio and extracted into a txt file for QC to be performed in R using the following scripts:

1. 01_bkgdcor_QC_prep.R - extracts the beta-values, methylation signals, and detection p-values from the GenomeStudio output

2. 02_bkgdcor_QC_CpGassoc.R - Samples with probe detection call rates <90% and those with an average intensity value of either <50% of the experiment-wide sample mean or <2,000 arbitrary units (AU) are excluded. Data points with probe detection p-values >0.001 are set to missing, and CpG sites with missing data for >10% of samples are excluded from analysis.  

3. 03_bkgdcor_beta_BMIQ_bySample.R - Probes that cross hybridize between autosomes and sex chromosomes (Chen 2013) are removed and Beta Mixture Quantile Normalization (BMIQ) is run.
