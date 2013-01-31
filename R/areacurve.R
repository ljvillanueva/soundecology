#
require("tuneR")
require("seewave")
require("pracma")

soundarea<-function(soundfile, min_freq=100, max_freq=10000){
	AUC_left=""
	AUC_right=""

	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		
		data_l <- meanspec(left,plot=FALSE)
		
		rows_1000 = dim(data_l)[1] / (left@samp.rate/2)
		min_row = min_freq * rows_1000
		max_row = max_freq * rows_1000
		
		data_l_segment <- data_l[min_row:max_row,]
		
		AUC_l = trapz(data_l_segment[,2])
		
		data_r <- meanspec(right,plot=FALSE)
		
		rows_1000 = dim(data_r)[1] / (right@samp.rate/2)
		min_row = min_freq * rows_1000
		max_row = max_freq * rows_1000
		
		data_r_segment <- data_r[min_row:max_row,]
		
		AUC_r = trapz(data_r_segment[,2])
		
		cat(paste(" AUC for ", min_freq, "-", max_freq, "Hz:\n", sep=""))
		
		cat(paste("  Left channel: ", AUC_l, "\n", sep=""))
		cat(paste("  Right channel: ", AUC_r, "\n\n", sep=""))
		rm(left, right)
	} else 
	{
		cat("\n This is a mono file.\n")
		
		#matrix of values
		left<-soundfile
		
		data_l <- meanspec(left,plot=FALSE)
		
		rows_1000 = dim(data_l)[1] / (left@samp.rate/2)
		min_row = min_freq * rows_1000
		max_row = max_freq * rows_1000
		
		data_l_segment <- data_l[min_row:max_row,]
		
		AUC_l = trapz(data_l_segment[,2])
		
		cat(paste(" AUC for ", min_freq, "-", max_freq, "Hz:\n", sep=""))
		
		cat(paste("  Mono channel: ", AUC_l, "\n", sep=""))
		rm(left)
	}
	invisible(list(AUC_left=AUC_l, AUC_right=AUC_r))
}
