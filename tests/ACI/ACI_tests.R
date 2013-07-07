#
#Test between ACI values from this package and the 
#  SoundscapeMeter v 1.0.14.05.2012 software, courtesy of A. Farina.
#

library(soundecology)

#testsound2_AciTot.txt <- sum the values for ACItot
#testsound2_AciIfTot.txt <- sum the values for ACIif

vals <- data.frame(matrix(ncol=3, nrow=250))

names(vals) <- c("filename", "ACI_orig", "ACI_soundecology")

files <- dir("tests/sounds/")

for (i in 1:length(files)){
#for (i in 1:10){
	cat(files[i])
	cat("\n")
	
	res_dir <- paste("tests/ACI/results/", strsplit(files[i], ".wav"), "/", sep="")
	res_file_ACItot <- paste(res_dir, strsplit(files[i], ".wav"), "_AciTot.txt", sep="")
	#res_file_ACIif <- paste(res_dir, strsplit(files[1], ".wav"), "_AciIfTot.txt", sep="")
	
	
	ACItot_vals <- read.csv(res_file_ACItot, header=FALSE, sep="\t", nrows=1)
	ACItot <- sum(ACItot_vals[1:length(ACItot_vals)-1])
	
	#ACIif_vals <- read.csv(res_file_ACIif, header=FALSE, sep="\t", nrows=1)
	#ACIif <- sum(ACIif_vals[1:length(ACIif_vals)-1])
	
	this_file <- paste("tests/sounds/", unlist(files[i]), sep="")
	soundfile <- readWave(this_file)
	ind_vals <- acoustic_complexity(soundfile, max_freq=22050, j=5, fft_w=512)
	
	vals[i, 1] <- files[i]
	vals[i, 2] <- ACItot
	vals[i, 3] <- ind_vals$AciTotAll_left
	}

cor.test(vals[,2], vals[,3])
wilcox.test(vals[,3], vals[,2], paired=TRUE)
