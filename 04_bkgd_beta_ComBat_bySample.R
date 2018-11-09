####################################################################################
# Background Corrected: ComBat 
####################################################################################

# Runs ComBat including PTSD and Gender (if necessary) as covariates

rm(list=ls())
library(sva)
library(impute)

####################################################################################
# Step 1: Load and Prep Data
####################################################################################

samps<-read.csv("SAMPLE SHEET.csv",row.names=1, stringsAsFactors=F)
load("bkgd_beta_BMIQ_bySample.Rdata")
beta.norm<-bmiq
rm(bmiq, design.v)
rownames(samps)<-samps$SentrixBarcode_Position
missing<-which(is.na(beta.norm)) # Creates flag for missing observations

# Note: if there are any beta values that are 0 or 1, they will have to be modified 
# before log transforming in order to prevent "-Inf"

beta.norm[which(beta.norm==1)]<-0.9999999
beta.norm[which(beta.norm==0)]<-0.0000001

beta.norm<-log2(beta.norm/(1-beta.norm)) # M-values
samps<-samps[colnames(beta.norm),]
all(rownames(samps)==colnames(beta.norm))

####################################################################################
# Step 2: Impute Data
####################################################################################

beta.imputed<-impute.knn(as.matrix(beta.norm))
beta.imputed<-beta.imputed$data
save(beta.imputed, file="bkgd_betaImputed_bySample_preComBat.Rdata")

####################################################################################
# Step 3: Run ComBat
####################################################################################

samps$SentrixBarcode_A<-as.factor(samps$SentrixBarcode_A)
samps$SentrixPosition_A<-as.factor(samps$SentrixPosition_A)

# model includes gender if necessary
mod<-model.matrix(~as.factor(PTSD)+as.factor(Gender), data=samps) 

# Run ComBat for Chip Effects, PTSD is a moderator
combat_beta<-ComBat(dat=beta.imputed,batch=samps$SentrixBarcode_A,mod=mod) 
save(combat_beta, file="bkgd_beta_bySample_ComBAT_chipAdj.Rdata")

# Run ComBat for Position Effects, PTSD is a moderator
combat_beta<-ComBat(dat=combat_beta,batch=samps$SentrixPosition_A, mod=mod)
save(combat_beta, file="bkgd_beta_bySample_ComBat_chipPosAdj.Rdata")

# Convert back to beta values
reversbeta<-2^combat_beta/(2^combat_beta + 1)
reversbeta[missing]<-NA # putting missing values back in
reversbeta<-round(reversbeta,digits=4)

save(reversbeta, missing, file="bkgd_beta_postComBat_bySample.Rdata")



