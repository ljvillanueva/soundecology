#Multiple sounds
#
# Script to process all .wav files in a directory and save
# the requested index to a .csv file.
#

multiple_sounds <- function(directory, resultfile, soundindex = c("ndsi", "acoustic_complexity", "acoustic_diversity", "bioacoustic_index", "H"), no_cores=1, ...){

	if (file.access(directory) == -1) {
		stop(paste("The directory specified does not exist or this user is not autorized to read it:\n    ", directory))
		}
	
	#is it windows?
	#if (.Platform$OS.type == "windows"){
	#	no_cores = 1
	#	}
	
	thismachine_cores <- detectCores()
	
	if (no_cores == 0){
		stop("Number of cores can not be 0.")
	}else if (no_cores < -1){
		stop("Number of cores can not be negative.")
	}else if (no_cores == "max"){
		no_cores = thismachine_cores
	}else if (no_cores == -1){
		no_cores = thismachine_cores - 1
	}else if (no_cores > thismachine_cores){
		#Don't try to use more than the number of cores in the machine
		warning(paste(" The number of cores to use can not be larger than the cores in the computer:", detectCores()), immediate.=TRUE)
		no_cores <- thismachine_cores
		}
		
	wav_files <- dir(path = directory, pattern="wav$", ignore.case=TRUE)
	if (length(wav_files)==0) {
		stop(paste("Could not find any .wav files in the specified directory:\n    ", directory))
	}
	

  
	if (soundindex == "bioacoustic_index"){
    
	  fileheader <- c("FILENAME,INDEX,FFT_W,MIN_FREQ,MAX_FREQ,LEFT_CHANNEL,RIGHT_CHANNEL\n")
    
		getindex <- function(soundfile, ...){
			library(soundecology)
			#Get args
			args <- list(...)
	
			if(!is.null(args[['min_freq']])) {
				min_freq = args[['min_freq']]
			}else{
				min_freq = formals(bioacoustic_index)$min_freq
			}
			if(!is.null(args[['max_freq']])) {
				max_freq = args[['max_freq']]
			}else{
				max_freq = formals(bioacoustic_index)$max_freq
			}
			if(!is.null(args[['fft_w']])) {
				fft_w = args[['fft_w']]
			}else{
				fft_w = formals(bioacoustic_index)$fft_w
			}
		
# 			if (file.access(resultfile) == -1) {
# 				cat("FILENAME,INDEX,FFT_W,MIN_FREQ,MAX_FREQ,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
# 			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			
			return_list <- bioacoustic_index(this_soundfile, ...)
			
			cat(paste(soundfile, ",", soundindex, ",", fft_w, ",", min_freq, ",", max_freq, ",", return_list$left_area, ",", return_list$right_area, "\n", sep=""), file=resultfile, append=TRUE)
			}
	}else if (soundindex == "acoustic_diversity"){
    
	  fileheader <- c("FILENAME,INDEX,MAX_FREQ,DB_THRESHOLD,FREQ_STEPS,LEFT_CHANNEL,RIGHT_CHANNEL\n")
	      
		getindex <- function(soundfile, ...){
			library(soundecology)
				
			#Get args
			args <- list(...)
			
			if(!is.null(args[['db_threshold']])) {
				db_threshold = args[['db_threshold']]
			}else{
				db_threshold = formals(acoustic_diversity)$db_threshold
			}
			if(!is.null(args[['max_freq']])) {
				max_freq = args[['max_freq']]
			}else{
				max_freq = formals(acoustic_diversity)$max_freq
			}
			if(!is.null(args[['freq_step']])) {
				freq_step = args[['freq_step']]
			}else{
				freq_step = formals(acoustic_diversity)$freq_step
			}
		
# 			if (file.access(resultfile) == -1) {
# 				cat("FILENAME,INDEX,MAX_FREQ,DB_THRESHOLD,FREQ_STEPS,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
# 			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_diversity(this_soundfile, ...)
		
			cat(paste(soundfile, ",", soundindex, ",", max_freq, ",", db_threshold, ",", freq_step, ",", return_list$shannon_left, ",", return_list$shannon_right, "\n", sep=""), file=resultfile, append=TRUE)
			}
	}else if (soundindex == "acoustic_complexity"){
    
	  fileheader <- c("FILENAME,INDEX,FFT_W,MAX_FREQ,J,LEFT_CHANNEL,RIGHT_CHANNEL\n")
    
		getindex <- function(soundfile, ...){
			library(soundecology)
			
			#Get args
			args <- list(...)
			
			if(!is.null(args[['max_freq']])) {
				max_freq = args[['max_freq']]
			}else{
				max_freq = formals(acoustic_complexity)$max_freq
			}
			if(!is.null(args[['j']])) {
				j = args[['j']]
			}else{
				j = formals(acoustic_complexity)$j
			}
			if(!is.null(args[['fft_w']])) {
				fft_w = args[['fft_w']]
			}else{
				fft_w = formals(acoustic_complexity)$fft_w
			}

# 			if (file.access(resultfile) == -1) {
# 				cat("FILENAME,INDEX,FFT_W,MAX_FREQ,J,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
# 			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- acoustic_complexity(this_soundfile, ...)
			
			cat(paste(soundfile, ",", soundindex, ",", fft_w, ",", max_freq, ",", j, ",", return_list$AciTotAll_left, ",", return_list$AciTotAll_right, "\n", sep=""), file=resultfile, append=TRUE)
			}
	}else if (soundindex == "ndsi"){
    
	  fileheader <- c("FILENAME,INDEX,FFT_W,ANTHRO_MIN,ANTHRO_MAX,BIO_MIN,BIO_MAX,HZ_INTERVAL,LEFT_CHANNEL,RIGHT_CHANNEL\n")
    
		getindex <- function(soundfile, ...){
			library(soundecology)
			
			#Get args
			args <- list(...)
			
			if(!is.null(args[['fft_w']])) {
				fft_w = args[['fft_w']]
			}else{
				fft_w = formals(ndsi)$fft_w
			}
			if(!is.null(args[['anthro_min']])) {
				anthro_min = args[['anthro_min']]
			}else{
				anthro_min = formals(ndsi)$anthro_min
			}
			if(!is.null(args[['anthro_max']])) {
				anthro_max = args[['anthro_max']]
			}else{
				anthro_max = formals(ndsi)$anthro_max
			}
			if(!is.null(args[['bio_min']])) {
				bio_min = args[['bio_min']]
			}else{
				bio_min = formals(ndsi)$bio_min
			}
			if(!is.null(args[['bio_max']])) {
				bio_max = args[['bio_max']]
			}else{
				bio_max = formals(ndsi)$bio_max
			}
			if(!is.null(args[['hz_interval']])) {
				hz_interval = args[['hz_interval']]
			}else{
				hz_interval = formals(ndsi)$hz_interval
			}
						
# 			if (file.access(resultfile) == -1) {
# 				cat("FILENAME,INDEX,FFT_W,ANTHRO_MIN,ANTHRO_MAX,BIO_MIN,BIO_MAX,HZ_INTERVAL,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
# 			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			return_list <- ndsi(this_soundfile, ...)
			cat(paste(soundfile, ",", soundindex, ",", fft_w, ",", anthro_min, ",", anthro_max, ",", bio_min, ",", bio_max, ",", hz_interval, ",", return_list$ndsi_left, ",", return_list$ndsi_right, "\n", sep=""), file=resultfile, append=TRUE)
			}
	}else if (soundindex == "H"){
    
	  fileheader <- c("FILENAME,INDEX,WL,ENVT,MSMOOTH,KSMOOTH,LEFT_CHANNEL,RIGHT_CHANNEL\n")
    
		getindex <- function(soundfile, ...){
			library(soundecology)
			
			#Get args
			args <- list(...)
			
			if(!is.null(args[['wl']])) {
				wl = args[['wl']]
			}else{
				wl = formals(H)$wl
			}
			if(!is.null(args[['envt']])) {
				envt = args[['envt']]
			}else{
				envt = formals(H)$envt
			}
			if(!is.null(args[['msmooth']])) {
				msmooth = args[['msmooth']]
			}else{
				msmooth = formals(H)$msmooth
			}
			if(!is.null(args[['ksmooth']])) {
				ksmooth = args[['ksmooth']]
			}else{
				ksmooth = formals(H)$ksmooth
			}
			
			library(soundecology)
# 			if (file.access(resultfile) == -1) {
# 				cat("FILENAME,INDEX,WL,ENVT,MSMOOTH,KSMOOTH,LEFT_CHANNEL,RIGHT_CHANNEL\n", file=resultfile, append=TRUE)
# 			}
			this_soundfile <- readWave(paste(directory, soundfile, sep=""))
			if (this_soundfile@stereo==TRUE) {
				left<-channel(this_soundfile, which = c("left"))
				right<-channel(this_soundfile, which = c("right"))
				left_res <- H(left, ...)
				right_res <- H(right, ...)
			}else{
				left<-channel(this_soundfile, which = c("left"))
				left_res <- H(left, ...)
				right_res <- NA
			}
			
			cat(paste(soundfile, ",", soundindex, ",", wl, ",", envt, ",", msmooth, ",", ksmooth, ",", left_res, ",", right_res, "\n", sep=""), file=resultfile, append=TRUE)
			}
		}

  
#open results file
sink(resultfile)
cat(fileheader)
#Done writing results  
sink()
  
#Use parallel?
if (no_cores>1){
	require(parallel)
	cat(paste("Running on ", length(wav_files), " files using ", no_cores, " cores", "\n\n", sep=""))
	
  cl <- makeCluster(no_cores, type = "PSOCK")
	
	res <- parLapply(cl, wav_files, getindex, ...)
	
	#pause to allow all to end
	Sys.sleep(1)
	
	stopCluster(cl)
}else{
  
	for (soundfile in wav_files){
		getindex(soundfile, ...)
		}
	}

	
}