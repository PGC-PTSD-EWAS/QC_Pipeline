###############################################################################################
# Minfi 450k Cell Count Estimation
###############################################################################################

library(minfi)
library(FlowSorted.Blood.450k)

baseDir =""

list.files(baseDir)
targets <- read.450k.sheet(baseDir)

rgSet<- read.450k.exp(base = baseDir, targets = targets)

cell.est<-estimateCellCounts(rgSet = rSet, compositeCellType = "Blood",
                             cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
                             returnAll = FALSE, meanPlot = FALSE, verbose = TRUE)

save.image("cell_estimates.Rdata")




