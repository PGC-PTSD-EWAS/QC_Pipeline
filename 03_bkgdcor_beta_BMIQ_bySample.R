####################################################################################
# Background Corrected: Removing Cross-hybridized probes and running BMIQ
####################################################################################

rm(list=ls())
library(wateRmelon)

####################################################################################
# Step 1: Loading beta matrix, annotation, and CpG sites to remove
####################################################################################

# 450K annotation for probe design
load("fullannotInd.rda")
workspace<-ls()
annot<-fullannot
rm(list=workspace, "workspace")

# CpGassoc QC'ed beta matrix
load("bkgdcor_beta_QC.Rdata")

# Cross-hybridized probes from Chen 2013
chen<-read.table("48639-non-specific-probes-Illumina450k.txt", header=T, sep="\t")

####################################################################################
# Step 2: Removing SNP-associated and Cross-hybridized probes
####################################################################################

# Creating list of probes to remove
probesRmv<-as.character(chen[, "TargetID"]) # list of probes to remove
probesRmv<-probesRmv[!is.na(match(probesRmv, rownames(beta.qc)))]
beta.qc<-beta.qc[-match(probesRmv, rownames(beta.qc)),]
rm(chen, probesRmv)

####################################################################################
# Step 3: Creating design.v. for BMIQ normalization.
####################################################################################

annot<-annot[rownames(beta.qc),]
beta.qc[which(beta.qc==0)]<-0.000001
beta.qc[which(beta.qc==1)]<-0.999999

design.v<-annot[, "INFINIUM_DESIGN_TYPE"]
design.v[design.v=="I"]<-1
design.v[design.v=="II"]<-2
design.v<-as.numeric(design.v)

####################################################################################
# Step 4: Running BMIQ
####################################################################################

# Checking order of beta and design matrices before running BMIQ
for(ii in 1:ncol(beta.qc)){
  beta.v<-beta.qc[, ii]
  sID<-colnames(beta.qc)[ii]
  bmiqTemp<-BMIQ(beta.v=beta.v, design.v=design.v, sampleID=sID) 
  beta.qc[,ii]<-bmiqTemp$nbeta
}

bmiq<-beta.qc
 
# Save list of all 8 outputs from BMIQ function
save(bmiq, design.v, file="bkgd_beta_BMIQ_bySample.Rdata")

rm(list=ls()[-match(c("bmiq", "design.v"), ls())])
