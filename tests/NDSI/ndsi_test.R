#
#Test between NDSI values from this package and the 
#  data from a Matlab script by Stuart Gage.
#

library(soundecology)

vals <- read.csv("calculated_index.csv")

cor.test(vals$Soundecology_NDSI, vals$MAT.AUCC_NDSI)

wilcox.test(vals$Soundecology_NDSI, vals$MAT.AUCC_NDSI, paired=TRUE)
plot(((vals$Soundecology_NDSI-vals$MAT.AUCC_NDSI)/vals$MAT.AUCC_NDSI)*100)
summary((vals$Soundecology_NDSI-vals$MAT.AUCC_NDSI)/vals$MAT.AUCC_NDSI)
boxplot((vals$Soundecology_NDSI-vals$MAT.AUCC_NDSI)/vals$MAT.AUCC_NDSI)



wilcox.test(vals$Soundecology_NDSI, vals$TLSP201306.0600.R_NDSI, paired=TRUE)
plot(((vals$Soundecology_NDSI-vals$TLSP201306.0600.R_NDSI)/vals$TLSP201306.0600.R_NDSI)*100)
summary((vals$Soundecology_NDSI-vals$TLSP201306.0600.R_NDSI)/vals$TLSP201306.0600.R_NDSI)
boxplot((vals$Soundecology_NDSI-vals$TLSP201306.0600.R_NDSI)/vals$TLSP201306.0600.R_NDSI)
