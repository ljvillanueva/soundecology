#
#Gage's Soundscape Index
# From REAL - Remote Environmental Assessment Laboratory
# http://www.real.msu.edu/
#
# Also from: Krause, Bernie, Stuart H. Gage, and Wooyeong Joo. 2011. 
#  Measuring and interpreting the temporal variability in the soundscape
#  at four places in Sequoia National Park. Landscape Ecol. DOI 10.1007/s10980-011-9639-6
#
# Also from: Kasten, Eric P., Stuart H. Gage, Jordan Fox, and Wooyeong Joo. 2012.
#  The Remote Environmental Assessment Laboratory’s Acoustic Library: An Archive for 
#  Studying Soundscape Ecology. Ecological Informatics 12: 50–67. doi:10.1016/j.ecoinf.2012.08.001.
#
# "The samples were converted to 22 kHz monaural. A normalized Power Spectrum
#  Density value (PSD in Watt/Hz) (Welch 1967) was
#  computed for each 1 kHz interval for all recordings
#  by running a script programmed by Gage and Joo
#  using MATLAB (Gilat 2004)."
#
# "Normalized Difference Soundscape Index (NDSI) is a numeric indicator of a
#  soundscape's relative biological composition to human disturbance based
#  on the amount of acoustic energy in different frequency bands. The distribution 
#  of sound frequencies in acoustic samples was computed to determine what 
#  types of sounds occurred (e.g., mechanical, biological or physical) at 
#  different frequencies. Normalized ratios of frequency levels were then 
#  calculated to evaluate the relative amounts of Biophony (biological sounds) 
#  or Anthrophony (human-induced sounds) in a set of samples. 
#  The equation for NDSI is: 
#  
#    NDSI = (β - α) / (β + α)
#  
#  where α, β represent the proportion of acoustic energy of Anthrophony and 
#  Biophony, respectively. The value of the index ranges from -1 to 1. If the 
#  value is negative, mechanical sounds dominate the soundscape. If the value 
#  approaches 1, most of the soundscape consists of biological sounds at the site."
#
#  "The NDSI is computed as the ratio of the sound intensity found in the 
#  frequencies where biological sounds (biophony) are most prevalent (2-8 kHz)
#  to the frequencies where mechanical sounds (technophony) are most prevalent
#  (1-2 kHz). NDSI has values in the range +1 to -1, with +1 indicating that 
#  a signal contains only biophony. As shown in the figure, NDSI has a larger 
#  biophonic component between 2100 and 0730 hours."
#
#
# REQUIRES the packages tuneR, seewave, pracma
#

ndsi <- function(soundfile, fft_w=512, anthro_min=1000, anthro_max=2000, bio_min=2000, bio_max=8000, hz_interval=1000){
	
	#Some general values
	#Get sampling rate
	samplingrate <- soundfile@samp.rate
	duration <- length(soundfile@left)/soundfile@samp.rate
	
	#Get Nyquist frequency in Hz
	nyquist_freq <- (samplingrate/2)
	
	#Check errors
	#if (max(anthro_max, bio_max)>nyquist_freq) {
	#	cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
	#	break
	#}
	
	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Calculating PSD... Please wait... \n")
		specA_left <- spec(left, f=samplingrate, wl=fft_w, plot=FALSE, PSD=TRUE)
		specA_right <- spec(right, f=samplingrate, wl=fft_w, plot=FALSE, PSD=TRUE)
		
		rm(left,right)
		
		#LEFT CHANNEL
		specA_rows <- dim(specA_left)[1]
		
		freq_per_row <- floor(specA_rows/nyquist_freq)
		
		delta_fl <- nyquist_freq / specA_rows

		anthro_vals_range <- anthro_max - anthro_min
		bio_vals_range <- bio_max - bio_min
		
		antro_rows <- anthro_vals_range / freq_per_row
		bio_rows <- bio_vals_range / freq_per_row
		
		#anthorphony_vals <- data.frame(matrix(NA, nrow = antro_rows, ncol = 2))
		#biophony_vals <- data.frame(matrix(NA, nrow = bio_rows, ncol = 2))
		
		anthro_min_row <- round(anthro_min / freq_per_row)
		anthro_max_row <- round(anthro_max / freq_per_row)
		bio_min_row <- round(bio_min / freq_per_row)
		bio_max_row <- round(bio_max / freq_per_row)
		
		#Left channel
		
		anthrophony_vals <- specA_left[anthro_min_row:anthro_max_row,1:2]
		biophony_vals <- specA_left[bio_min_row:bio_max_row,1:2]
		
		anthrophony_auc <- trapz(anthrophony_vals[,1], anthrophony_vals[,2])
		biophony_auc <- trapz(biophony_vals[,1], biophony_vals[,2])
		
		NDSI_left = (biophony_auc - anthrophony_auc) / (biophony_auc + anthrophony_auc)

		#Right channel
		anthrophony_vals <- specA_right[anthro_min_row:anthro_max_row,1:2]
		biophony_vals <- specA_right[bio_min_row:bio_max_row,1:2]
		
		anthrophony_auc <- trapz(anthrophony_vals[,1], anthrophony_vals[,2])
		biophony_auc <- trapz(biophony_vals[,1], biophony_vals[,2])
		
		NDSI_right = (biophony_auc - anthrophony_auc) / (biophony_auc + anthrophony_auc)
				
		cat("\n NDSI left channel: ")
		cat(NDSI_left)
		cat("\n")
		cat(" NDSI right channel: ")
		cat(NDSI_right)
		cat("\n\n")
	} else 
	{
		cat("\n This is a mono file.\n")
		left<-channel(soundfile, which = c("left"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Calculating PSD... Please wait... \n")
		specA_left <- spec(left, f=samplingrate, wl=fft_w, plot=FALSE, PSD=TRUE)
		
		rm(left)
		
		#LEFT CHANNEL
		specA_rows <- dim(specA_left)[1]
		
		freq_per_row <- floor(specA_rows/nyquist_freq)
		
		anthro_vals_range <- anthro_max - anthro_min
		bio_vals_range <- bio_max - bio_min
		
		antro_rows <- anthro_vals_range / freq_per_row
		bio_rows <- bio_vals_range / freq_per_row
		
		#anthorphony_vals <- data.frame(matrix(NA, nrow = antro_rows, ncol = 2))
		#biophony_vals <- data.frame(matrix(NA, nrow = bio_rows, ncol = 2))
		
		anthro_min_row <- round(anthro_min / freq_per_row)
		anthro_max_row <- round(anthro_max / freq_per_row)
		bio_min_row <- round(bio_min / freq_per_row)
		bio_max_row <- round(bio_max / freq_per_row)
		
		#Left channel
		
		anthrophony_vals <- specA_left[anthro_min_row:anthro_max_row,1:2]
		biophony_vals <- specA_left[bio_min_row:bio_max_row,1:2]
		
		anthrophony_auc <- trapz(anthrophony_vals[,1], anthrophony_vals[,2])
		biophony_auc <- trapz(biophony_vals[,1], biophony_vals[,2])
		
		NDSI_left = (biophony_auc - anthrophony_auc) / (biophony_auc + anthrophony_auc)
		
		#Right channel
		NDSI_right = NA
		
		cat("\n NDSI: ")
		cat(NDSI_left)
		cat("\n\n")
	}
	invisible(list(ndsi_left=NDSI_left, ndsi_right=NDSI_right))
}
