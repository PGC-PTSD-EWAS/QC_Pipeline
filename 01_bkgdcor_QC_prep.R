
####################################################################################
# Background Corrected Data Prep for QC
####################################################################################

# Input file contains the background corrected beta value, methylation signals
# and detection p-values for all CpG sites

####################################################################################
# Loading Data
####################################################################################

tab <- read.table("GENOME STUDIO OUTPUT.txt", header = TRUE, 
	nrows = 5, sep="\t", check.names=F)

classes <- sapply(tab, class)
classes[which(classes=="integer")]<-"numeric"	
rm(tab)

tab <- read.table("GENOME STUDIO OUTPUT.txt", header = TRUE, 
	sep="\t", check.names=F, colClasses = classes)
rownames(tab)<-tab[, "TargetID"]
rm(classes)

####################################################################################
# Step 1: Extract Beta matrix
####################################################################################

beta<-tab[, grep("AVG_Beta", colnames(tab))]
str(beta)
dim(beta)
betam<-as.matrix(beta)
head(betam)[,1:4]
head(beta)[,1:4]
all(beta==betam, na.rm=T)
all(which(is.na(beta))==which(is.na(betam)))
beta<-betam
colnames(beta)<-gsub(".AVG_Beta", "", colnames(beta))
str(beta)

save(beta, file="bkgdcor_beta_RAW.Rdata")

rm(betam, beta)

####################################################################################
# Step 2: Extract Signal A
####################################################################################

sigA<-tab[, grep("Signal_A", colnames(tab))]
str(sigA)
dim(sigA)
sigAm<-as.matrix(sigA)
head(sigAm)[,1:4]
head(sigA)[,1:4]
all(sigA==sigAm, na.rm=T)
all(which(is.na(sigA))==which(is.na(sigAm)))
sigA<-sigAm
colnames(sigA)<-gsub(".Signal_A", "", colnames(sigA))
str(sigA)

save(sigA, file="bkgdcor_signalA_RAW.Rdata")

rm(sigA, sigAm)

####################################################################################
# Step 3: Extract Signal B
####################################################################################

sigB<-tab[, grep("Signal_B", colnames(tab))]
str(sigB)
dim(sigB)
sigBm<-as.matrix(sigB)
head(sigBm)[,1:4]
head(sigB)[,1:4]
all(sigB==sigBm, na.rm=T)
all(which(is.na(sigBm))==which(is.na(sigB)))
sigB<-sigBm
colnames(sigB)<-gsub(".Signal_B", "", colnames(sigB))
str(sigB)

save(sigB, file="bkgdcor_signalB_RAW.Rdata")

rm(sigB, sigBm)

####################################################################################
# Step 4: Extract Detection p-values
####################################################################################

detP<-tab[, grep("Detection", colnames(tab))]
str(detP)
dim(detP)
detPm<-as.matrix(detP)
head(detP)[,1:4]
head(detPm)[,1:4]
all(detP==detPm, na.rm=T)
all(which(is.na(detP))==which(is.na(detPm)))
detP<-detPm
colnames(detP)<-gsub(".Detection Pval", "", colnames(detP))

save(detP, file="bkgdcor_detectP_RAW.Rdata")
