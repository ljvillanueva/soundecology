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

getSplitFile <- function(time,file,j,len) {

  if(time+j>len){
    return ()
  }
  else{
    return (readWave(file,from=time,to=time+j,units="seconds"))
  }
  #return readWave(file,from = start, to=end, units = "seconds")
}

acousticComplexitySplit <- function
(
  soundfile,
  fileLength,
  minFreq = NA,
  maxFreq = NA,
  j = 5,
  fft_w = 512,
  split = NA,
  matrix = TRUE,
  bands = TRUE
) {

  #file length given in milliseconds
  lengthSeconds = fileLength / 1000
  #cutting length to size
  lengthSeconds = lengthSeconds - lengthSeconds%%j
  timeVector = seq(0,lengthSeconds,j)

  splitSound = lapply(timeVector,getSplitFile,file=soundfile,j=j,len=lengthSeconds)
  splitSound[sapply(splitSound, is.null)] <- NULL
  results = sapply(splitSound,acoustic_complexity_new,minFreq,maxFreq,j=j,fft_w=fft_w,matrix=matrix,bands=bands)

  numRows = 5
  if(matrix){
    numRows = numRows + 2
  }
  if(bands){
    numRows = numRows + 2
  }
  numCols = dim(results)[2]
  print(numCols)
  numRows = dim(results)[1]
  print(numRows)
  
  return(
    list(
      aciTotAllLeft <- sum(unlist(results[1,1:numCols])),
      aciTotAllRight <- sum(unlist(results[2,1:numCols])),
      aciTotLeftByMin <- round( (aciTotAllLeft / lengthSeconds) * 60, 2),
      aciTotRightByMin <- round( (aciTotAllRight / lengthSeconds) * 60, 2),
      aciOverTimeLeft <- c(unlist(results[5,1:numCols])),
      aciOverTimeRight <- c(unlist(results[6,1:numCols]))
    )
  )
  }



acoustic_complexity_new <- function(
  soundfile,
  minFreq = NA,
  maxFreq = NA,
  j = 5,
  fft_w = 512,
  matrix = TRUE,
  bands = TRUE) {
  # Helper Functions
  check_param <- function(soundfile, minFreq, maxFreq, j, fft_w) {
    # Temp Values
    MIN <- minFreq
    MAX <- maxFreq
    J <- j
    FFT_W <- fft_w

    # test arguments
    if (is.na(MAX)) {
      MAX <- soundfile@samp.rate / 2
      cat(paste("\n maxFreq not set, using value of:", MAX, "\n\n"))
    } else if (MAX > nyquistFreq) {
      cat(
        paste(
          "\n WARNING:
          The maximum acoustic frequency that this file can use is ",
           nyquistFreq,
           "Hz. But the script was set to measure up to ",
           MAX,
           "Hz. The value of maxFreq was changed to ",
           nyquistFreq, ".\n\n",
           sep = ""
         )
      )
      MAX <- nyquistFreq
      # break
    }

    if (is.na(MIN)) {
      MIN <- 0
      cat(paste("\n minFreq not set, using value of:", MIN, "\n\n"))
    }

    if (is.numeric(as.numeric(MIN))) {
      MIN <- as.numeric(MIN)
    } else {
      stop(" minFreq is not a number.")
    }

    if (is.numeric(as.numeric(MAX))) {
      MAX <- as.numeric(MAX)
    } else {
      stop(" maxFreq is not a number.")
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
          round( duration %% j,digits=5),
          " seconds have been discarded from output.",
          sep = ""
        )
      )
    }
    if (is.numeric(as.numeric(FFT_W))) {
      FFT_W <- as.numeric(FFT_W)
    } else {
      stop(" fft_w is not a number.")
    }

    eval.parent(substitute(minFreq <- MIN))
    eval.parent(substitute(maxFreq <- MAX))
    eval.parent(substitute(j <- J))
    eval.parent(substitute(fft_w <- FFT_W))
  }

  # function that gets the difference of values
  getD <- function(index, spectrum, freqRow, minCol, maxCol, iPerJ) {
    minCol <- index * iPerJ - iPerJ + 1
    maxCol <- index * iPerJ

    return(
      sum(
        abs(
          spectrum[
                    freqRow,
                    minCol:(maxCol - 1)
                  ] -
          spectrum[
                    freqRow,
                    (minCol + 1):maxCol
                  ]
          )
       )
    )
  }

  # Vectorized function to get sum of amps within range/freq_band
  sumV <- function(index, spectrum, qIndex) {
    minCol <- index * iPerJ - iPerJ + 1
    maxCol <- index * iPerJ
    return(sum(spectrum[qIndex, minCol:maxCol]))
  }


  # Setting Duration,`` Window length, sampling rate, nyquist freq respectively
  duration <- length(soundfile@left) / soundfile@samp.rate
  wlen <- fft_w
  samplingrate <- soundfile@samp.rate
  nyquistFreq <- (samplingrate / 2)

  # Checking validity of parameters
  check_param(soundfile, minFreq, maxFreq, j, fft_w)



  # Stereo file
  if (soundfile@stereo == TRUE) {
    cat("\n This is a stereo file. Results will be given for each channel.\n")
    # TODO: Make sure the new code is the same as this code
    left <- soundfile@left
    right <- soundfile@right

    # matrix of values
    cat("\n Calculating index. Please wait... \n\n")
    specLeft <- spectro(
                      left,
                      f = samplingrate,
                      wl = wlen,
                      plot = FALSE,
                      norm = TRUE,
                      dB = NULL,
                      scale = FALSE,
                      wn = "hamming"
                  )
    specAmpLeft <- specLeft$amp
    specLeftFreq <- specLeft$freq
    rm(specLeft)

    minFreqByThousand <- minFreq / 1000
    maxFreqByThousand <- maxFreq / 1000

    minFreqIndex <- which(
                        abs(specLeftFreq - minFreqByThousand) ==
                        min(abs(specLeftFreq - minFreqByThousand))
                      )
    maxFreqIndex <- which(
                        abs(specLeftFreq - maxFreqByThousand) ==
                        min(abs(specLeftFreq - maxFreqByThousand))
                      )

    if (minFreqIndex < 1) {
      minFreqIndex <- 1
    }

    if (maxFreqIndex > dim(specAmpLeft)[1]) {
      maxFreqIndex <- dim(specAmpLeft)[1] - 1
    }


    numSpecRows <- dim(specAmpLeft)[1]
    numSpecCols <- dim(specAmpLeft)[2]

    deltaTK <- (length(soundfile@left) / soundfile@samp.rate) / numSpecCols
    rm(left)
    numJ <- floor(duration / j)

    # Number of values, in each row, for each j period (no. of columns)
    iPerJ <- floor(j / deltaTK)

    aciLeftVals <- rep(numJ)
    aciComplexityLeft <- rep(0, numJ)

    aciRightVals <- rep(numJ)
    aciComplexityRight <- rep(0, numJ)

    if (matrix) {
      aciLeftMatrix <- matrix(nrow = numSpecRows, ncol = numJ)
      aciRightMatrix <- matrix(nrow = numSpecRows, ncol = numJ)
    }
    if (bands) {
      aciFlLeftVector <- rep(numJ)
      aciFlRightVector <- rep(numJ)
    } else {
      totalLeft <- 0
      totalRight <- 0
    }

    # create index for vectorized function
    indexes <- 1:numJ

    # Left channel
    # For each frequency bin fl
    for (qIndex in 1:numSpecRows) {

      # For each j period of time
      D <- sapply(
              indexes,
              getD,
              spectrum = specAmpLeft,
              freqRow = qIndex,
              minCol = minCol,
              maxCol = maxCol,
              iPerJ = iPerJ
            )
      sumI <- sapply(
                  indexes,
                  sumV,
                  spectrum = specAmpLeft,
                  qIndex = qIndex
                )

      aciLeftVals <- D / sumI
      aciComplexityLeft <- aciComplexityLeft + aciLeftVals

      if (matrix) {
        aciLeftMatrix[qIndex, ] <- aciLeftVals
      }

      if (bands) {
        aciFlLeftVector[qIndex] <- sum(aciLeftVals)
      } else {
        totalLeft <- totalLeft + sum(aciLeftVals)
      }

    }

    if (bands) {
      aciTotLeft <- sum(aciFlLeftVector)
    } else {
      aciTotLeft <- totalLeft
    }
    # Clean Up Amp matrix
    rm(specAmpLeft)


    # Retrieve spectrogram of right side
    specRight <- spectro(
                      right,
                      f = samplingrate,
                      wl = wlen,
                      plot = FALSE,
                      norm = TRUE,
                      dB = NULL,
                      scale = FALSE,
                      wn = "hamming"
                    )

    specAmpRight <- specRight$amp[minFreqIndex:maxFreqIndex, ]

    rm(specRight)
    rm(right)
    # Right channel
    # For each frequency bin fl
    for (qIndex in 1:numSpecRows) {
      # For each j period of time

      D <- sapply(
              indexes,
              getD,
              spectrum = specAmpRight,
              freqRow = qIndex,
              minCol = minCol,
              maxCol = maxCol,
              iPerJ = iPerJ
            )
      sumI <- sapply(
                  indexes,
                  sumV,
                  spectrum = specAmpRight,
                  qIndex = qIndex
                )
      aciRightVals <- D / sumI
      aciComplexityRight <- aciComplexityRight + aciRightVals
      if (matrix) {
        aciRightMatrix[qIndex, ] <- aciRightVals
      }
      if (bands) {
        aciFlRightVector[qIndex] <- sum(aciRightVals)
      } else {
        totalRight <- totalRight + sum(aciRightVals)
      }
    }

    if (bands) {
      aciTotRight <- sum(aciFlRightVector)
    } else {
      aciTotRight <- totalRight
    }

    aciTotLeftByMin <- round( (aciTotLeft / duration) * 60, 2)
    aciTotRightByMin <- round( (aciTotRight / duration) * 60, 2)

    cat(
      paste(
        "Acoustic Complexity Index (total):\n",
        "   Left channel: ",
        sep = ""
      )
    )
    cat(aciTotLeft)
    cat(paste("\n", "   Right channel: ", sep = ""))
    cat(aciTotRight)
    cat("\n\n")
    if (duration > 60) {
      cat(
        paste(
          "  Acoustic Complexity Index (by minute):\n",
          "   Left channel: ",
          sep = ""
        )
      )
      cat(aciTotLeftByMin)
      cat(paste("\n", "   Right channel: ", sep = ""))
      cat(aciTotRightByMin)
      cat("\n\n")
    }
  } else {
    cat("\n This is a mono file.\n")

    left <- channel(soundfile, which = c("left"))

    # matrix of values
    cat("\n Calculating index. Please wait... \n\n")
    specLeft <- spectro(
                    left,
                    f = samplingrate,
                    wl = wlen,
                    plot = FALSE,
                    norm = TRUE,
                    dB = NULL,
                    scale = FALSE,
                    wn = "hamming"
                )

    specAmpLeft <- specLeft$amp
    specLeftFreq <- specLeft$freq
    rm(specLeft)

    minFreqByThousand <- minFreq / 1000
    maxFreqByThousand <- maxFreq / 1000

    minFreqIndex <- which(
                        abs(specLeftFreq - minFreqByThousand) ==
                        min(abs(specLeftFreq - minFreqByThousand))
                      )
    maxFreqIndex <- which(
                        abs(specLeftFreq - maxFreqByThousand) ==
                        min(abs(specLeftFreq - maxFreqByThousand))
                      )

    specAmpLeft <- specAmpLeft[minFreqIndex:maxFreqIndex, ]
    rm(left)

    # LEFT CHANNEL
    numSpecRows <- dim(specAmpLeft)[1]
    numSpecCols <- dim(specAmpLeft)[2]
    #
    # 		freq_per_row <- numSpecRows/nyquistFreq
    #
    # 		max_row <- round(max_freq * freq_per_row)
    #
    # 		specAmpLeft <- specAmpLeft[1:max_row,]
    # 		numSpecRows <- dim(specAmpLeft)[1]

    deltaTK <- (length(soundfile@left) / soundfile@samp.rate) / numSpecCols

    numJ <- floor(duration / j)
    # q <- numSpecRows
    # m <- floor(duration / j)

    # Number of values, in each row, for each j period (no. of columns)
    iPerJ <- floor(j / deltaTK)

    aciLeftVals <- rep(NA, numJ)
    aciRightVals <- rep(NA, numJ)

    if (matrix) {
      aciLeftMatrix <- matrix(nrow = numSpecRows, ncol = numJ)
      aciRightMatrix <- matrix(nrow = numSpecRows, ncol = numJ)
    }
    if (bands) {
      aciFlLeftVector <- rep(numJ)
      aciFlRightVector <- rep(numJ)
    } else {
      totalLeft <- 0
      totalRight <- 0
    }

    #create index vector
    indexes <- 1:numJ
    # Left channel
    # For each frequency bin fl
    for (qIndex in 1:numSpecRows) {


      D <- sapply(
              indexes,
              getD,
              spectrum = specAmpLeft,
              freqRow = qIndex,
              minCol = minCol,
              maxCol = maxCol,
              iPerJ = iPerJ
            )
      sumI <- sapply(
                  indexes,
                  sumV,
                  spectrum = specAmpLeft,
                  qIndex = qIndex
                )
      aciLeftVals <- D / sumI

      aciComplexityLeft <- aciComplexityLeft + aciLeftVals
      if (matrix) {
        aciLeftMatrix[qIndex, ] <- aciLeftVals
      }
      if (bands) {
        aciFlLeftVector[qIndex] <- sum(aciLeftVals)
      } else {
        totalLeft <- totalLeft + sum(aciLeftVals)
      }
    }

    if (bands) {
      aciTotLeft <- sum(aciFlLeftVector)
    } else {
      aciTotLeft <- totalLeft
    }

    aciTotLeftByMin <- round( (aciTotLeft / duration) * 60, 2)

    aciTotRight <- NA
    aciTotRightByMin <- NA

    cat("  Acoustic Complexity Index (total): ")
    cat(aciTotLeft)
    cat("\n\n")
    if (duration > 60) {
      cat("  Acoustic Complexity Index (by minute): ")
      cat(aciTotLeftByMin)
      cat("\n\n")
    }
  }

  return <- list(
    aciTotAllL = aciTotLeft,
    aciTotAllR = aciTotRight,
    aciTotAllByMinL = aciTotLeftByMin,
    aciTotAllByMinR = aciTotRightByMin,
    aciOverTimeL = aciComplexityLeft,
    aciOverTimeR = aciComplexityRight
  )

  if (matrix) {
    matricies <- list(
      aciLeftMatrix = aciLeftMatrix,
      aciRightMatrix = aciRightMatrix
    )

    return <- c(return, matricies)
  }

  if (bands) {
    bands <- list(
      aciFlLeftVals = aciFlLeftVector,
      aciFlRightVals = aciFlRightVector
    )
    return <- c(return, bands)
  }

  invisible(return)
}
