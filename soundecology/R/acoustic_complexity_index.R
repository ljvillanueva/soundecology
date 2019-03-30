#
# Acoustic Complexity Index
# From: N. Pieretti, A. Farina, D. Morri. 2011. A new methodology to infer
# Contributer: Ot Gabaldon Torrents, Josh Pollmann
# the singing activity of an avian community: The Acoustic Complexity Index (ACI).
# Ecological Indicators 11: 868-873.
#
# Tested with SoundscapeMeter 1.0.14.05.2012, courtesy of A. Farina
#
# TODO: Clean up code, SPlit files on J, look over other indicies

get_split_time <- function(j) {
  z <- 0
  while (j * z < 5) {
    z <- z + 1
  }
  return(z * j)
}

acoustic_complexity_split <- function
(
  soundfile,
  min_freq = NA,
  max_freq = NA,
  j = 5, fft_w = 512,
  split = NA,
  matrix = TRUE,
  bands = TRUE
) {
  duration <- length(soundfile@left) / soundfile@samp.rate
  if (duration > 10) {
    if (is.na(split)) {
      split_time <- get_split_time(j)
    }

    num_splits <- floor(duration / split_time)
    split_times <- 1:num_splits * split_time
    split_times <- c(split_times, duration)

    print(split_times)
  } else {
    return(
      acoustic_complexity_new(
        soundfile,
        min_freq,
        max_freq,
        j,
        fft_w,
        matrix,
        bands)
    )
  }
}


acoustic_complexity_new <- function(
  soundfile,
  min_freq = NA,
  max_freq = NA,
  j = 5,
  fft_w = 512,
  matrix = TRUE,
  bands = TRUE) {
  # Helper Functions
  check_param <- function(soundfile, min_freq, max_freq, j, fft_w) {
    # Temp Values
    MIN <- min_freq
    MAX <- max_freq
    J <- j
    FFT_W <- fft_w

    # test arguments
    if (is.na(MAX)) {
      MAX <- soundfile@samp.rate / 2
      cat(paste("\n max_freq not set, using value of:", MAX, "\n\n"))
    } else if (MAX > nyquist_freq) {
      cat(
        paste(
          "\n WARNING:
          The maximum acoustic frequency that this file can use is ",
           nyquist_freq,
           "Hz. But the script was set to measure up to ",
           MAX,
           "Hz. The value of max_freq was changed to ",
           nyquist_freq, ".\n\n",
           sep = ""
         )
      )
      MAX <- nyquist_freq
      # break
    }

    if (is.na(MIN)) {
      MIN <- 0
      cat(paste("\n min_freq not set, using value of:", MIN, "\n\n"))
    }

    if (is.numeric(as.numeric(MIN))) {
      MIN <- as.numeric(MIN)
    } else {
      stop(" min_freq is not a number.")
    }

    if (is.numeric(as.numeric(MAX))) {
      MAX <- as.numeric(MAX)
    } else {
      stop(" max_freq is not a number.")
    }

    if (is.numeric(as.numeric(J))) {
      J <- as.numeric(J)
    } else {
      stop(" j is not a number.")
    }
    # TODO: This is temporary fix,
    #       make sure to elongate ACI value acccording to cut value.
    if (duration < J) {
      J <- duration
      cat(paste(
        "J value, ", J, ", exceeds file duration, ", duration, ".",
        sep = ""
      ))
    }
    if (duration %% J != 0) {
      cat(
        paste(
          ( duration %% j),
          " seconds have been discarded from output",
          sep = ""
        )
      )
    }
    if (is.numeric(as.numeric(FFT_W))) {
      FFT_W <- as.numeric(FFT_W)
    } else {
      stop(" fft_w is not a number.")
    }

    eval.parent(substitute(min_freq <- MIN))
    eval.parent(substitute(max_freq <- MAX))
    eval.parent(substitute(j <- J))
    eval.parent(substitute(fft_w <- FFT_W))
  }

  # function that gets the difference of values
  get_d <- function(index, spectrum, freq_row, min_col, max_col, I_per_j) {
    min_col <- index * I_per_j - I_per_j + 1
    max_col <- index * I_per_j

    return(
      sum(
        abs(
          spectrum[
                    freq_row,
                    min_col:(max_col - 1)
                  ] -
          spectrum[
                    freq_row,
                    (min_col + 1):max_col
                  ]
          )
       )
    )
  }

  # Vectorized function to get sum of amps within range/freq_band
  sum_v <- function(index, spectrum, q_index) {
    min_col <- index * I_per_j - I_per_j + 1
    max_col <- index * I_per_j
    return(sum(spectrum[q_index, min_col:max_col]))
  }


  # Setting Duration,`` Window length, sampling rate, nyquist freq respectively
  duration <- length(soundfile@left) / soundfile@samp.rate
  wlen <- fft_w
  samplingrate <- soundfile@samp.rate
  nyquist_freq <- (samplingrate / 2)
  # Checking validity of parameters
  check_param(soundfile, min_freq, max_freq, j, fft_w)



  # Stereo file
  if (soundfile@stereo == TRUE) {
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    # TODO: Make sure the new code is the same as this code
    left <- soundfile@left
    right <- soundfile@right

    # matrix of values
    cat("\n Calculating index. Please wait... \n\n")
    spec_left <- spectro(
                      left,
                      f = samplingrate,
                      wl = wlen,
                      plot = FALSE,
                      norm = TRUE,
                      dB = NULL,
                      scale = FALSE,
                      wn = "hamming"
                  )

    spec_amp_left <- spec_left$amp
    spec_left_freq <- spec_left$freq
    rm(spec_left)

    min_freq1k <- min_freq / 1000
    max_freq1k <- max_freq / 1000

    which_min_freq <- which(
                        abs(spec_left_freq - min_freq1k) ==
                        min(abs(spec_left_freq - min_freq1k))
                      )
    which_max_freq <- which(
                        abs(spec_left_freq - max_freq1k) ==
                        min(abs(spec_left_freq - max_freq1k))
                      )

    if (which_min_freq < 1) {
      which_min_freq <- 1
    }

    if (which_max_freq > dim(spec_amp_left)[1]) {
      which_max_freq <- dim(spec_amp_left)[1] - 1
    }


    spec_amp_rows <- dim(spec_amp_left)[1]
    spec_amp_cols <- dim(spec_amp_left)[2]

    delta_tk <- (length(soundfile@left) / soundfile@samp.rate) / spec_amp_cols

    no_j <- floor(duration / j)

    # Number of values, in each row, for each j period (no. of columns)
    I_per_j <- floor(j / delta_tk)

    ACI_left_vals <- rep(no_j)
    ACI_complex_left_vector <- rep(0, no_j)

    ACI_right_vals <- rep(no_j)
    ACI_complex_right_vector <- rep(0, no_j)

    if (matrix) {
      ACI_left_matrix <- matrix(nrow = spec_amp_rows, ncol = no_j)
      ACI_right_matrix <- matrix(nrow = spec_amp_rows, ncol = no_j)
    }
    if (bands) {
      ACI_fl_left_vector <- rep(no_j)
      ACI_fl_right_vector <- rep(no_j)
    } else {
      total_left <- 0
      total_right <- 0
    }

    # create index for vectorized function
    index_v <- 1:no_j

    # Left channel
    # For each frequency bin fl
    for (q_index in 1:spec_amp_rows) {

      # For each j period of time
      D <- sapply(
              index_v,
              get_d,
              spectrum = spec_amp_left,
              freq_row = q_index,
              min_col = min_col,
              max_col = max_col,
              I_per_j = I_per_j
            )
      sum_I <- sapply(
                  index_v,
                  sum_v,
                  spectrum = spec_amp_left,
                  q_index = q_index
                )
      ACI_left_vals <- D / sum_I

      ACI_complex_left_vector <- ACI_complex_left_vector + ACI_left_vals
      if (matrix) {
        ACI_left_matrix[q_index, ] <- ACI_left_vals
      }
      if (bands) {
        ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
      } else {
        total_left <- total_left + sum(ACI_left_vals)
      }
    }

    if (bands) {
      ACI_tot_left <- sum(ACI_fl_left_vector)
    } else {
      ACI_tot_left <- total_left
    }
    # Clean Up Amp matrix
    rm(spec_amp_left)


    # Retrieve spectrogram of right side
    spec_right <- spectro(
                      right,
                      f = samplingrate,
                      wl = wlen,
                      plot = FALSE,
                      norm = TRUE,
                      dB = NULL,
                      scale = FALSE,
                      wn = "hamming"
                    )

    spec_amp_right <- spec_right$amp[which_min_freq:which_max_freq, ]

    rm(spec_right)
    rm(right)
    # Right channel
    # For each frequency bin fl
    for (q_index in 1:spec_amp_rows) {
      # For each j period of time

      D <- sapply(
              index_v,
              get_d,
              spectrum = spec_amp_right,
              freq_row = q_index,
              min_col = min_col,
              max_col = max_col,
              I_per_j = I_per_j
            )
      sum_I <- sapply(
                  index_v,
                  sum_v,
                  spectrum = spec_amp_right,
                  q_index = q_index
                )
      ACI_right_vals <- D / sum_I
      ACI_complex_right_vector <- ACI_complex_right_vector + ACI_right_vals
      if (matrix) {
        ACI_right_matrix[q_index, ] <- ACI_right_vals
      }
      if (bands) {
        ACI_fl_right_vector[q_index] <- sum(ACI_right_vals)
      } else {
        total_right <- total_right + sum(ACI_right_vals)
      }
    }

    if (bands) {
      ACI_tot_right <- sum(ACI_fl_right_vector)
    } else {
      ACI_tot_right <- total_right
    }

    ACI_tot_left_by_min <- round( (ACI_tot_left / duration) * 60, 2)
    ACI_tot_right_by_min <- round( (ACI_tot_right / duration) * 60, 2)

    cat(
      paste(
        "Acoustic Complexity Index (total):\n",
        "   Left channel: ",
        sep = ""
      )
    )
    cat(ACI_tot_left)
    cat(paste("\n", "   Right channel: ", sep = ""))
    cat(ACI_tot_right)
    cat("\n\n")
    if (duration > 60) {
      cat(
        paste(
          "  Acoustic Complexity Index (by minute):\n",
          "   Left channel: ",
          sep = ""
        )
      )
      cat(ACI_tot_left_by_min)
      cat(paste("\n", "   Right channel: ", sep = ""))
      cat(ACI_tot_right_by_min)
      cat("\n\n")
    }
  } else {
    cat("\n This is a mono file.\n")

    left <- channel(soundfile, which = c("left"))

    # matrix of values
    cat("\n Calculating index. Please wait... \n\n")
    spec_left <- spectro(
                    left,
                    f = samplingrate,
                    wl = wlen,
                    plot = FALSE,
                    norm = TRUE,
                    dB = NULL,
                    scale = FALSE,
                    wn = "hamming"
                )

    spec_amp_left <- spec_left$amp

    min_freq1k <- min_freq / 1000
    max_freq1k <- max_freq / 1000

    which_min_freq <- which(
                        abs(spec_left$freq - min_freq1k) ==
                        min(abs(spec_left$freq - min_freq1k))
                      )
    which_max_freq <- which(
                        abs(spec_left$freq - max_freq1k) ==
                        min(abs(spec_left$freq - max_freq1k))
                      )

    spec_amp_left <- spec_amp_left[which_min_freq:which_max_freq, ]
    rm(spec_left)

    rm(left)

    # LEFT CHANNEL
    spec_amp_rows <- dim(spec_amp_left)[1]
    spec_amp_cols <- dim(spec_amp_left)[2]
    #
    # 		freq_per_row <- spec_amp_rows/nyquist_freq
    #
    # 		max_row <- round(max_freq * freq_per_row)
    #
    # 		spec_amp_left <- spec_amp_left[1:max_row,]
    # 		spec_amp_rows <- dim(spec_amp_left)[1]

    fl <- rep(NA, spec_amp_rows)
    delta_fl <- (max_freq - min_freq) / spec_amp_rows
    delta_tk <- (length(soundfile@left) / soundfile@samp.rate) / spec_amp_cols

    no_j <- floor(duration / j)
    # q <- spec_amp_rows
    # m <- floor(duration / j)

    # Number of values, in each row, for each j period (no. of columns)
    I_per_j <- floor(j / delta_tk)

    ACI_left_vals <- rep(NA, no_j)
    ACI_fl_left_vector <- rep(NA, no_j)
    ACI_left_matrix <- data.frame(matrix(NA, nrow = spec_amp_rows, ncol = no_j))

    ACI_right_vals <- rep(NA, no_j)
    ACI_fl_right_vector <- rep(NA, no_j)
    ACI_right_matrix <- data.frame(matrix(NA, nrow = spec_amp_rows, ncol = no_j))

    # Left channel
    # For each frequency bin fl
    for (q_index in 1:spec_amp_rows) {

      # For each j period of time
      for (j_index in 1:no_j) {
        min_col <- j_index * I_per_j - I_per_j + 1
        max_col <- j_index * I_per_j

        D <- get_d(spec_amp_left, q_index, min_col, max_col)
        sum_I <- sum(spec_amp_left[q_index, min_col:max_col])
        ACI_left_vals[j_index] <- D / sum_I
        ACI_left_matrix[q_index, j_index] <- D / sum_I
      }

      ACI_fl_left_vector[q_index] <- sum(ACI_left_vals)
    }

    ACI_tot_left <- sum(ACI_fl_left_vector)
    ACI_tot_left_by_min <- round( (ACI_tot_left / duration) * 60, 2)

    ACI_tot_right <- NA
    ACI_tot_right_by_min <- NA

    cat("  Acoustic Complexity Index (total): ")
    cat(ACI_tot_left)
    cat("\n\n")
    if (duration > 60) {
      cat("  Acoustic Complexity Index (by minute): ")
      cat(ACI_tot_left_by_min)
      cat("\n\n")
    }
  }

  return <- list(
    AciTotAll_left = ACI_tot_left,
    AciTotAll_right = ACI_tot_right,
    AciTotAll_left_bymin = ACI_tot_left_by_min,
    AciTotAll_right_bymin = ACI_tot_right_by_min,
    ACI_complex_right_vector = ACI_complex_right_vector,
    ACI_complex_left_vector = ACI_complex_left_vector
  )

  if (matrix) {
    matricies <- list(
      aci_left_matrix = ACI_left_matrix,
      aci_right_matrix = ACI_right_matrix
    )

    return <- c(return, matricies)
  }

  if (bands) {
    bands <- list(
      aci_fl_left_vals = ACI_fl_left_vector,
      aci_fl_right_vals = ACI_fl_right_vector
    )
    return <- c(return, bands)
  }

  invisible(return)
}
