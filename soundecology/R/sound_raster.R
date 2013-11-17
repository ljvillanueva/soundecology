#Acoustic Diversity Index from Villanueva-Rivera \emph{et al.} 2011. 
# The ADI is calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default -50 dBFS). The ADI is the result of the Shannon index applied to these bins.

sound_raster <- function(wavfile=NA, max_freq=10000, wav_directory=FALSE){

  if (is.na(wavfile)==TRUE && wav_directory==FALSE){
    stop(" You have to provide a filename, in the argument wavfile, or a directory as the argument wav_directory.\n\n")
  }

  max_freq <- as.numeric(max_freq)

  soundfile<-readWave(wavfile)
  
  #Get Nyquist frequency in Hz
  nyquist_freq<-(soundfile@samp.rate/2)
  
  if (max_freq > nyquist_freq) {
    stop(paste("\n WARNING: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
  }
  
  #window length for the spectro and spec functions
  #to keep each row every 10Hz
  #Frequencies and seconds covered by each
  freq_per_row = 10
  wlen = soundfile@samp.rate/freq_per_row
  
	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n This is a stereo file. A raster will be written for each channel. Please wait...\n")
		left<-channel(soundfile, which = c("left"))
		right<-channel(soundfile, which = c("right"))
		rm(soundfile)
		
		#matrix of values
		specA_left <- spectro(left, f=samplingrate, wl=wlen, plot=FALSE)$amp
		specA_right <- spectro(right, f=samplingrate, wl=wlen, plot=FALSE)$amp
		
		rm(left, right)
		
		#left channel
		max_freq_row <- max_freq/10
		sound_map_left <- specA_left[1:max_freq_row,]
		
		nrows <- dim(sound_map_left)[1]
		ncols <- dim(sound_map_left)[2]
		
		#Create container array for db values
		sound_map <- array(dim=dim(sound_map_left))
		
		#Using k to write the new array with lower freqs in the bottom
		k=dim(sound_map)[1]
		for (i in seq(from=1, to=dim(sound_map)[1], by=1)){
		  for (j in seq(from=1, to=dim(sound_map)[2], by=1)){
		    sound_map[k,j]=round(sound_map_left[i,j], digits = 4)
		  }
		  k=k-1
		}
		
    #ascii filename
		ascii_raster_left <- paste(strsplit(wavfile,".wav"), "_wav_left.asc", sep="")
		#Write header
		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep=""), file=ascii_raster_left, append=FALSE)
		
		#Write data
		write.table(sound_map, file=ascii_raster_left, append=TRUE, row.names=FALSE, col.names=FALSE, sep=" ")
		
	
    
    #Right channel
		sound_map_right <- specA_right[1:max_freq_row,]
				
		#Create container array for db values
		sound_map <- array(dim=dim(sound_map_right))
		
		#Using k to write the new array with lower freqs in the bottom
		k <- dim(sound_map)[1]
		for (i in seq(from=1, to=dim(sound_map)[1], by=1)){
		  for (j in seq(from=1, to=dim(sound_map)[2], by=1)){
		    sound_map[k,j] <- round(sound_map_right[i,j], digits = 4)
		  }
		  k <- k-1
		}
		
		#ascii filename
		ascii_raster_right <- paste(strsplit(wavfile,".wav"), "_wav_right.asc", sep="")
		#Write header
		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep=""), file = ascii_raster_right, append=FALSE)
		
		#Write data
		write.table(sound_map, file = ascii_raster_right, append = TRUE, row.names = FALSE, col.names = FALSE, sep = " ")
		
		cat(paste("\n  ASCII raster files: \n    ", ascii_raster_left, "\n    ", ascii_raster_right, "\n\n", sep = ""))
    
	} else {

		cat("\n This is a mono file. A single raster will be written for this file. Please wait...\n")
		left<-channel(soundfile, which = c("left"))
		rm(soundfile)
		
		#matrix of values
		specA_left <- spectro(left, f=samplingrate, wl=wlen, plot=FALSE)$amp
		
		rm(left)
		
		#left channel
		max_freq_row <- max_freq/10
		sound_map_left <- specA_left[1:max_freq_row,]
		
		nrows <- dim(sound_map_left)[1]
		ncols <- dim(sound_map_left)[2]
		
		#Create container array for db values
		sound_map <- array(dim=dim(sound_map_left))
		
		#Using k to write the new array with lower freqs in the bottom
		k=dim(sound_map)[1]
		for (i in seq(from=1, to=dim(sound_map)[1], by=1)){
		  for (j in seq(from=1, to=dim(sound_map)[2], by=1)){
		    sound_map[k,j]=round(sound_map_left[i,j], digits = 4)
		  }
		  k=k-1
		}
		
		#ascii filename
		ascii_raster_left <- paste(strsplit(wavfile,".wav"), "_wav_left.asc", sep="")
		#Write header
		cat(paste("ncols         ", ncols, "\nnrows         ", nrows, "\nxllcorner     0.0\nyllcorner     0.0\ncellsize      1\nNODATA_value  -9999\n", sep=""), file=ascii_raster_left, append=FALSE)
		
		#Write data
		write.table(sound_map, file=ascii_raster_left, append=TRUE, row.names=FALSE, col.names=FALSE, sep=" ")
		
    cat(paste("\n  ASCII raster file: \n    ", ascii_raster_left, "\n\n", sep=""))
    
	}
}
