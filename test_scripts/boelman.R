#Bioacoustic index
#
#From the paper:
# Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion
#   resistance in Hawaii: bioacoustics, field surveys, and airborne
#   remote sensing. Ecol Applications 17(8):2137-44.
#
# Inspired on Matlab code provided by NT Boelman. Still need to check if
#  results are comparable.
# Boelman et al. 2007 used min_freq=2000, max_freq=8000, fft_w=1024

bioacoustic_index<-function(soundfile, min_freq=2000, max_freq=8000, fft_w=1024){
	#Get sampling rate
	samplingrate<-soundfile@samp.rate
	freq_per_row = 10
	wlen=samplingrate/freq_per_row
	
	#Get Nyquist frequency in Hz
	nyquist_freq<-(samplingrate/2)
	
	if (max_freq>nyquist_freq) {
		cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
		break
		}
	
	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n Stereo file.\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")
		specA_left <- meanspec(left, f=samplingrate, wl=fft_w, plot=FALSE, dB="max0")[,2]
		rm(left)
		specA_right <- meanspec(right, f=samplingrate, wl=fft_w, plot=FALSE, dB="max0")[,2]
		rm(right)
		
		#Select rows
		rows_width = length(specA_left) / nyquist_freq
		
		min_row = min_freq * rows_width
		max_row = max_freq * rows_width
		
		#specA_left_segment <- specA_left[min_row:max_row,]
		specA_left_segment <- specA_left[min_row:max_row]
		specA_right_segment <- specA_right[min_row:max_row]
		
		freq_range <- max_freq - min_freq
		
		specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)
		left_area <- moredB(specA_left_segment_normalized) * freq_range
		
		specA_right_segment_normalized <- specA_right_segment - min(specA_right_segment)
		right_area <- moredB(specA_right_segment_normalized) * freq_range
		
		cat("\n")
		cat(paste(" Bioacoustic Index (Frequency range: ", min_freq, "-", max_freq, " Hz; FFT window of ", fft_w, "):\n", sep=""))
		
		cat("  Left channel: ")
		cat(left_area)
		cat("\n  Right channel: ")
		cat(right_area)
		cat("\n\n")
	} else 
	{
		cat("\n Mono file.\n")
		left<-channel(soundfile, which = c("left"))
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")
		#specA_left <- spectro(left, f=samplingrate, wl=1024, plot=FALSE)$amp
		specA_left <- meanspec(left, f=samplingrate, wl=fft_w, plot=FALSE, dB="max0")[,2]
		#cat(specA_left)
		rm(left)
		
		#left_avg <- rowMeans(specA_left) #This may be wrong, these are logs. May need to convert to amplitude to calculate
		#left_avg <- meandB(specA_left)
		
		#Select rows
		#rows_width = dim(specA_left)[1] / nyquist_freq
		rows_width = length(specA_left) / nyquist_freq
		
		min_row = min_freq * rows_width
		max_row = max_freq * rows_width
		
		#specA_left_segment <- specA_left[min_row:max_row,]
		specA_left_segment <- specA_left[min_row:max_row]
		
		freq_range <- max_freq - min_freq
		
		specA_left_segment_normalized <- specA_left_segment - min(specA_left_segment)
		
		left_area <- moredB(specA_left_segment_normalized) * freq_range
		
		cat("\n")
		cat(paste(" Bioacoustic Index (Frequency range: ", min_freq, "-", max_freq, " Hz; FFT window of ", fft_w, "):\n", sep=""))
		
		cat("  Mono channel: ")
		cat(left_area)
		AUC_r = trapz(specA_left_segment_2)
		cat("\n\n")
		right_area <- NA
	}
	invisible(list(left_area=left_area, right_area=right_area))
}

require("tuneR")
require("seewave")


soundfile <- 
