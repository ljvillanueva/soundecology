#Multiple sounds
#
# Script to process all .wav files in a directory and save
# the requested index to a .csv file.
#

multiple_sounds <- function(directory, resultfile, soundindex, ...){

	if (file.access(directory) == -1) {
		stop(paste("The directory specified does not exist or this user is not autorized to read it:\n    ", directory))
		}
	
	wav_files <- dir(path = directory, pattern="wav$", ignore.case=TRUE)
	if (length(wav_files)==0) {
		stop(paste("Could not find any .wav files in the specified directory:\n    ", directory))
	}
	
	for (soundfile in wav_files){
		if (soundindex == "bioacoustic_index"){
			if (file.access(resultfile) == -1) {
				cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- bioacoustic_index(this_soundfile, ...)
			cat(paste(soundfile, ",", soundindex, ",", return_list$left_area, ",", return_list$right_area, "\n", sep=""), file=resultfile, append=TRUE)
		}else if (soundindex == "acoustic_diversity"){
			if (file.access(resultfile) == -1) {
				cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_diversity(this_soundfile, ...)
			cat(paste(soundfile, ",", soundindex, ",", return_list$shannon_left, ",", return_list$shannon_right, "\n", sep=""), file=resultfile, append=TRUE)
		}else if (soundindex == "acoustic_complexity"){
			if (file.access(resultfile) == -1) {
				cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_complexity(this_soundfile, ...)
			cat(paste(soundfile, ",", soundindex, ",", return_list$AciTotAll_left, ",", return_list$AciTotAll_right, "\n", sep=""), file=resultfile, append=TRUE)
		}else if (soundindex == "ndsi"){
			if (file.access(resultfile) == -1) {
				cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- ndsi(this_soundfile, ...)
			cat(paste(soundfile, ",", soundindex, ",", return_list$ndsi_left, ",", return_list$ndsi_right, "\n", sep=""), file=resultfile, append=TRUE)
		}else if (soundindex == "H"){
			if (file.access(resultfile) == -1) {
				cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			if (this_soundfile@stereo==TRUE) {
				left<-channel(this_soundfile, which = c("left"))
				right<-channel(this_soundfile, which = c("right"))
				left_res <- H(left, ...)
				right_res <- H(right, ...)
				cat(paste(soundfile, ",", soundindex, ",", left_res, ",", right_res, "\n", sep=""), file=resultfile, append=TRUE)
			}else{
				left<-channel(this_soundfile, which = c("left"))
				left_res <- H(left, ...)
				right_res <- NA
				cat(paste(soundfile, ",", soundindex, ",", left_res, ",", right_res, "\n", sep=""), file=resultfile, append=TRUE)
			}
		}else{
			stop("  The value of soundindex was not recognized.")
		}
	}
}
