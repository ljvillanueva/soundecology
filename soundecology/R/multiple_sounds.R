#Multiple sounds
#
# Script to process all .wav files in a directory and save
# the requested index to a .csv file.
#

multiple_sounds <- function(directory, filename, soundindex, ...){

	if (file.access(directory) == -1) {
		stop(paste("The directory specified does not exist or this user is not autorized to read it:\n    ", directory))
		}
	
	wav_files <- dir(path = directory, pattern="wav$", ignore.case=TRUE)
	if (length(wav_files)==0) {
		stop(paste("Could not find any .wav files in the specified directory:\n    ", directory))
	}
	
	return_list <- NA
	
	if (file.access(filename) == -1) {
		cat("FILENAME,INDEX,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=filename, append=TRUE)
		}
	
	for (soundfile in wav_files){
		if (soundindex == "bioacoustic_index"){
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- bioacoustic_index(this_soundfile, ...)
			cat(paste(soundfile, soundindex, return_list$left_area, return_list$right_area, sep=","), file=filename, append=TRUE)
			cat("\n", file=filename, append=TRUE)
		}else if (soundindex == "acoustic_diversity"){
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_diversity(this_soundfile, ...)
			cat(paste(soundfile, soundindex, return_list$shannon_left, return_list$shannon_right, sep=","), file=filename, append=TRUE)
			cat("\n", file=filename, append=TRUE)
		}else if (soundindex == "acoustic_complexity"){
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_complexity(this_soundfile, ...)
			cat(paste(soundfile, soundindex, return_list$AciTotAll_left, return_list$AciTotAll_right, sep=","), file=filename, append=TRUE)
			cat("\n", file=filename, append=TRUE)
		}else if (soundindex == "ndsi"){
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- ndsi(this_soundfile, ...)
			cat(paste(soundfile, soundindex, return_list$ndsi_left, return_list$ndsi_right, sep=","), file=filename, append=TRUE)
			cat("\n", file=filename, append=TRUE)
		}
	}
}
