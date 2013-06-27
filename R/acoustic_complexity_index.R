#
#Acoustic Complexity Index
#From: N. Pieretti ∗ , A. Farina, D. Morri. 2011. A new methodology to infer
# the singing activity of an avian community: The Acoustic Complexity Index (ACI).
# Ecological Indicators 11: 868–873.
#
#Tested with SoundscapeMeter.1.0.14.05.2012, courtesy of A. Farina
#
acoustic_complexity <- function(soundfile, max_freq=22050, j=5, fft_w=512){
	
	#function that gets the difference of values
	# probably a very inefficient way, find better way
	get_d <- function(spectrum, freq_row, min_col, max_col){
		D = 0
		for (k in min_col:(max_col-1)) {
			D = D + abs(spectrum[freq_row,k] - spectrum[freq_row,k + 1])
			}
			 		
		return(D)
		}
	
	#Some general values
	#Get sampling rate
	samplingrate <- soundfile@samp.rate
	duration <- length(soundfile@left)/soundfile@samp.rate
	
	#Get Nyquist frequency in Hz
	nyquist_freq <- (samplingrate/2)
	if (max_freq>nyquist_freq) {
		cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", nyquist_freq, "Hz. But the script was set to measure up to ", max_freq, "Hz.\n\n", sep=""))
		break
	}
	
	#window length for the spectro and spec functions
	wlen = fft_w
	
	#Stereo file
	if (soundfile@stereo==TRUE) {
		cat("\n This is a stereo file. Results will be given for each channel.\n")
		left <- channel(soundfile, which = c("left"))
		right <- channel(soundfile, which = c("right"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")
		specA_left <- spectro(left, f=samplingrate, wl=wlen, plot=FALSE, norm=FALSE, dB=NULL, scale=FALSE)$amp
		specA_right <- spectro(right, f=samplingrate, wl=wlen, plot=FALSE, norm=FALSE, dB=NULL, scale=FALSE)$amp
		
		rm(left,right)
		
		#LEFT CHANNEL
		specA_rows <- dim(specA_left)[1]
		specA_cols <- dim(specA_left)[2]
		
		fl <- rep(NA, specA_rows)
		delta_fl <- nyquist_freq / specA_rows
		delta_tk <- (length(soundfile@left)/soundfile@samp.rate) / specA_cols
				
		m <- floor(duration / j)
		q <- specA_rows
		no_j <- floor(duration / j)
		
		#Number of values, in each row, for each j period (no. of columns)
		I_per_j <- floor(j/delta_tk)
		
		ACI_left_vals <- rep(NA, no_j)
		ACI_fl_left_vector <- rep(NA, m)
		ACI_left_matrix <- data.frame(matrix(NA, nrow = q, ncol = m))
		
		ACI_right_vals <- rep(NA, no_j)
		ACI_fl_right_vector <- rep(NA, m)
		ACI_right_matrix <- data.frame(matrix(NA, nrow = q, ncol = m))
		
		#Left channel
		#For each frequency bin fl
		for (q_index in 1:q) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j
								
				D <- get_d(specA_left, q_index, min_col, max_col)
				sum_I <- sum(specA_left[q_index,min_col:max_col])
				ACI_left_vals[j_index] <- D / sum_I
				ACI_left_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
			} 
		
		ACI_tot_left <- sum(ACI_fl_left_vector)
		
		#Right channel
		#For each frequency bin fl
		for (q_index in 1:q) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j				
				
				D <- get_d(specA_right, q_index, min_col, max_col)
				sum_I <- sum(specA_right[q_index,min_col:max_col])
				ACI_right_vals[j_index] <- D / sum_I
				ACI_right_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_right_vector[q_index] <- sum(ACI_right_vals)
			} 
		
		ACI_tot_right <- sum(ACI_fl_right_vector)
		
		cat(" ACI total left channel: ")
		cat(ACI_tot_left)
		cat("\n")
		cat(" ACI total right channel: ")
		cat(ACI_tot_right)
		cat("\n\n")
		
	} else 
	{
		cat("\n This is a mono file.\n")
		
		left<-channel(soundfile, which = c("left"))
		rm(soundfile)
		
		#matrix of values
		cat("\n Getting values from spectrogram... Please wait... \n")
		specA_left <- spectro(left, f=samplingrate, wl=wlen, plot=FALSE, norm=FALSE, dB=NULL, scale=FALSE)$amp
		
		rm(left)
		
		#LEFT CHANNEL
		specA_rows <- dim(specA_left)[1]
		specA_cols <- dim(specA_left)[2]
		
		fl <- rep(NA, specA_rows)
		delta_fl <- nyquist_freq / specA_rows
		delta_tk <- (length(soundfile@left)/soundfile@samp.rate) / specA_cols
		
		m <- floor(duration / j)
		q <- specA_rows
		no_j <- floor(duration / j)
		
		#Number of values, in each row, for each j period (no. of columns)
		I_per_j <- floor(j/delta_tk)
		
		ACI_left_vals <- rep(NA, no_j)
		ACI_fl_left_vector <- rep(NA, m)
		ACI_left_matrix <- data.frame(matrix(NA, nrow = q, ncol = m))
		
		ACI_right_vals <- rep(NA, no_j)
		ACI_fl_right_vector <- rep(NA, m)
		ACI_right_matrix <- data.frame(matrix(NA, nrow = q, ncol = m))
		
		#Left channel
		#For each frequency bin fl
		for (q_index in 1:q) {
			
			#For each j period of time
			for (j_index in 1:no_j) {
				min_col <- j_index * I_per_j - I_per_j + 1
				max_col <- j_index * I_per_j
								
				D <- get_d(specA_left, q_index, min_col, max_col)
				sum_I <- sum(specA_left[q_index,min_col:max_col])
				ACI_left_vals[j_index] <- D / sum_I
				ACI_left_matrix[q_index, j_index] <- D / sum_I
			}
			
			ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
		} 
		
		ACI_tot_left <- sum(ACI_fl_left_vector)
		
		ACI_tot_right <- NA
		
		cat(" ACI total: ")
		cat(ACI_tot_left)
		cat("\n\n")
	}
	
	invisible(list(AciTotAll_left=ACI_tot_left, AciTotAll_right=ACI_tot_right, 
				   aci_fl_left_vals=ACI_fl_left_vector, aci_fl_right_vals=ACI_fl_right_vector,
				   aci_left_matrix=ACI_left_matrix, aci_right_matrix=ACI_right_matrix))
}
