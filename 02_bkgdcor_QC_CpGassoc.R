####################################################################################
# Background Corrected Analysis CpGassoc QC
####################################################################################

# CpGassoc QC requires a matrix of beta values, the unmethylated (signal A) and
# methylated (signal B) signals, and the detection P-values from the GenomeStudio
# output. The script below separates out the values for each from the All samples
# tab file

library(CpGassoc)

load("bkgdcor_beta_RAW.Rdata")

load("bkgdcor_signalA_RAW.Rdata")

load("bkgdcor_detectP_RAW.Rdata")

load("bkgdcor_signalB_RAW.Rdata")

# Checking the order of the matrices. These should all be true
all(colnames(beta)==colnames(detP))
all(rownames(beta)==rownames(detP))
all(colnames(beta)==colnames(sigA))
all(rownames(beta)==rownames(sigA))
all(colnames(beta)==colnames(sigB))
all(rownames(beta)==rownames(sigB))

beta.qc<-cpg.qc(beta.orig = beta, siga = sigA, sigb = sigB, pval = detP, 
	p.cutoff = 0.001, cpg.miss=0.1)

save(beta.qc, file="bkgdcor_beta_QC.Rdata")






