#Boelman index
#Translated from her Matlab code
#From the paper:
# Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion
#   resistance in Hawaii: bioacoustics, field surveys, and airborne
#   remote sensing. Ecol Applications 17(8):2137-44.
#

boelman_index<-function(soundfile, min_freq=180, max_freq=8500, fft_w=1024){
	
	#Get sampling rate
	samplingrate<-soundfile@samp.rate
	
	#Get Nyquist frequency in Hz
	nyquist_freq<-(samplingrate/2)
	
	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")
		specA_left <- spectro(left, f=samplingrate, wl=wlen, plot=FALSE)$amp
		specA_right <- spectro(right, f=samplingrate, wl=wlen, plot=FALSE)$amp
		
		rm(left,right)
		
		if (max_freq>nyquist_freq) {
			cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
			break
		}
		
		Freq<-seq(from=0, to=max_freq-freq_step, by=freq_step)
		
		#LEFT CHANNEL
		
		#Score=seq(from=0, to=0, length=length(Freq))
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j]=getscore(specA_left, Freq[j], (Freq[j]+freq_step), db_threshold, freq_per_row)
		}
		
		left_vals=Score
		
		Score1=0
		for (i in 1:length(Freq)) {
			Score1=Score1 + (Score[i] * log(Score[i]+0.0000001))
		}
		
		#Average
		Score_left=(-(Score1))/length(Freq)
		
		
		#RIGHT CHANNEL
		
		#Score=seq(from=0, to=0, length=length(Freq))
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j]=getscore(specA_right, Freq[j], (Freq[j]+freq_step), db_threshold, freq_per_row)
		}
		
		right_vals=Score
		
		Score1=0
		for (i in 1:length(Freq)) {
			Score1=Score1 + (Score[i] * log(Score[i]+0.0000001))
		}
		
		#Average
		Score_right=(-(Score1))/length(Freq)
		
		
		cat(" ==============================================\n")
		cat(paste(" Results (with a dB threshold of ", db_threshold, ")\n\n", sep=""))
		
		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
		cat("Frequency range (Hz), left channel proportion, right channel proportion\n")
		for (j in seq(length(Freq),1,by=-1)) {
			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), ",", round(right_vals[j],6), "\n", sep=""))
		}
		
		cat("\n Plot of proportions in each band: \n\n")
		cat("  Left channel\n")
		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq),1,by=-1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j]+freq_step),"",sep="")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1,this_row_space,by=1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep="")
			}
			
			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
			temp_val=round(left_vals[j],2)*20
			if (temp_val>0){
				for (i in 1:temp_val) {
					cat("*")
				}
			}
			cat("\n")
			rm(temp_val)
		}
		
		cat("\n  Right channel\n")
		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq),1,by=-1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j]+freq_step),"",sep="")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1,this_row_space,by=1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep="")
			}
			
			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
			
			temp_val=round(right_vals[j],2)*20
			if (temp_val>0){
				for (i in 1:temp_val) {
					cat("*")
				}
			}
			cat("\n")
			rm(temp_val)
		}
		
		cat("\n  Acoustic diversity (Shannon's Index): \n")
		cat(paste("   Left channel: ", round(Score_left,6), "\n", sep=""))
		cat(paste("   Right channel: ", round(Score_right,6), "\n\n", sep=""))
		left_adi_return = round(Score_left,6)
		right_adi_return = round(Score_right,6)
		
		cat("  Band Eveness (Gini coefficient): \n")
		cat(paste("   Left channel: ", round(Gini(left_vals),6), "\n", sep=""))
		cat(paste("   Right channel: ", round(Gini(right_vals),6), "\n\n", sep=""))
		left_gini_return = round(Gini(left_vals),6)
		right_gini_return = round(Gini(right_vals),6)
		
	} else 
	{
		cat("\n This is a mono file.\n")
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")

		left<-channel(soundfile, which = c("left"))
		rm(soundfile)
		
		if (max_freq>nyquist_freq) {
			cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
			break
		}
		
		data_l <- meanspec(left, f=samplingrate, wl=fft_w, dB="max0", plot=FALSE)
		
		rows_1000 = dim(data_l)[1] / (left@samp.rate/2)
		min_row = min_freq * rows_1000
		max_row = max_freq * rows_1000
		
		data_l_segment <- data_l[min_row:max_row,]
		
		freq_range = max_freq - min_freq
		
		data_l_segment[,2]
		
		data_l_segment2 <- data_l_segment[,2] - min(data_l_segment[,2])
		
		soundarea <- sum(data_l_segment2 * freq_range)
		
		#Way too big number... should be around 20-40
		
		
		
		
		AUC_l = trapz(data_l_segment[,2])
		
		
		
		
		
		
		Freq<-seq(from=0, to=max_freq-freq_step, by=freq_step)
		
		#Score=seq(from=0, to=0, length=length(Freq))
		Score <- rep(NA, length(Freq))
		
		for (j in 1:length(Freq)) {
			Score[j]=getscore(specA_left, Freq[j], (Freq[j]+freq_step), db_threshold, freq_per_row)
		}
		
		left_vals=Score
		
		Score1=0
		for (i in 1:length(Freq)) {
			Score1=Score1 + (Score[i] * log(Score[i]+0.0000001))
		}
		
		#Average
		Score_left=(-(Score1))/length(Freq)
		
		
		cat(" ==============================================\n")
		cat(paste(" Results (with a dB threshold of ", db_threshold, ")\n\n", sep=""))
		
		cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
		cat("Frequency range (Hz), proportion\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq),1,by=-1)) {
			cat(paste(Freq[j], "-", (Freq[j]+freq_step), ",", round(left_vals[j],6), "\n", sep=""))
		}
		
		cat("\n Plot of proportions in each band: \n\n")
		cat("   Freq. range (Hz) |--------------------|\n")
		
		#printed in inverse order to keep the low frequencies in the bottom, like in a spectrogram
		for (j in seq(length(Freq),1,by=-1)) {
			this_row_name <- paste(Freq[j], "-", (Freq[j]+freq_step),"",sep="")
			this_row_size <- nchar(this_row_name)
			this_row_space <- 17 - this_row_size
			
			this_row_spaces = ""
			
			for (f in seq(1,this_row_space,by=1)) {
				this_row_spaces = paste(this_row_spaces, " ", sep="")
			}
			
			cat(paste("   ", this_row_name, this_row_spaces, "|", sep=""))
			temp_val=round(left_vals[j],2)*20
			if (temp_val>0){
				for (i in 1:temp_val) {
					cat("*")
				}
			}
			cat("\n")
			rm(temp_val)
		}
		
		cat("\n  Acoustic diversity (Shannon's Index): ")
		cat(paste(round(Score_left,6), "\n", sep=""))
		left_adi_return = round(Score_left,6)
		right_adi_return = 0
		
		cat("  Band Eveness (Gini coefficient): ")
		cat(paste(round(Gini(left_vals),6), "\n", sep=""))
		left_gini_return = round(Gini(left_vals),6)
		right_gini_return = 0
	}
	invisible(list(adi_left=left_adi_return, adi_right=right_adi_return, gini_left=left_gini_return, gini_right=right_gini_return))
}
