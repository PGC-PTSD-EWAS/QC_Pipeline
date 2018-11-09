####################################################################################
# Barfield Population Stratification Probes within 1 bp
####################################################################################

####################################################################################
# NOTES:
# This code uses the non-BMIQ normalized beta matrix that has gone through the 
# the CpGassoc QC removing probes with missing data >10%, removing samples with
# probe detection call rates <90% and those with an average intensity value of
# either <50% of the experiment-wide sample mean or <2,000 AUs.
####################################################################################

load("bkgdcor_beta_QC.Rdata")

####################################################################################
# Start Barfield Code:

### Sample code for use with location.based.pc() function

### This code can be used as an example or can be directly modified by the user to compute "location-based" 
### principal components as described in Barfield et al. (2014) Gen Epid, in press. 

### file.pathname should be a character string containing the local pathname for the annotation file
### Annotation files (based on 1000 Genomes Project data) should be downloaded from
### http://genetics.emory.edu/research/?assetID=2161

### Modify this line to indicate correct pathname
file.pathname = "cpgs_within_1bp_of_TGP_SNP.Rdata"

### beta.obj is the users data; should be a data.frame or matrix of beta values or M-values with: 
###### one row per CpG site, one column per sample
###### row names = CpG names, column names = sample names  

### Modify this line to incorporate user methylation data
beta.obj = beta.qc

### Load and call location.based.pc() function to compute principal components
### This function will select only CpG sites in the location-based annotation file,
### set missing values to the CpG-site-average, and compute principal components

source('location.based.pc.R')
pc <- location.based.pc(beta.obj,file.pathname)

save(pc, file="bkgdcor_beta_QC_1pc.Rdata")

### The returned result will be a princomp object.  The principal components will be available as pc$loadings

### Principal components can be incorporated as covariates in regression analysis
### Samples will be sorted in the same order as beta.obj, but if using with other data, make sure IDs match up!

### Example of using principal components in regression analysis

#top10pc <- pc$loadings[,1:10]
##library(CpGassoc)
#cpgassoc_result <- cpg.assoc(beta.val=beta.obj, indep=PHENOTYPE, cov=top10pc)
