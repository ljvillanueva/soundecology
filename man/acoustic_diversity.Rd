\name{acoustic_diversity}
\alias{acoustic_diversity}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Acoustic Diversity of a sound file
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
acoustic_diversity(soundfile, save_to_text = FALSE, max_freq = 10000, db_threshold = -50, freq_step = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{soundfile}{
%%     ~~Describe \code{soundfile} here~~
}
  \item{save_to_text}{
%%     ~~Describe \code{save_to_text} here~~
}
  \item{max_freq}{
%%     ~~Describe \code{max_freq} here~~
}
  \item{db_threshold}{
%%     ~~Describe \code{db_threshold} here~~
}
  \item{freq_step}{
%%     ~~Describe \code{freq_step} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (soundfile, save_to_text = FALSE, max_freq = 10000, 
    db_threshold = -50, freq_step = 1000) 
{
    getscore <- function(spectrum, minf, maxf, db, freq_row) {
        miny <- round((minf)/freq_row)
        maxy <- round((maxf)/freq_row)
        subA = spectrum[miny:maxy, ]
        index1 <- length(subA[subA > db])/length(subA)
        return(index1)
    }
    samplingrate <- soundfile@samp.rate
    nyquist_freq <- (samplingrate/2)
    freq_per_row = 10
    wlen = samplingrate/freq_per_row
    if (soundfile@stereo == TRUE) {
        cat("\n This is a stereo file. Results will be given for each channel.\n")
        left <- channel(soundfile, which = c("left"))
        right <- channel(soundfile, which = c("right"))
        rm(soundfile)
        cat("\n Getting values from spectrogram... Please wait... \n")
        specA_left <- spectro(left, f = samplingrate, wl = wlen, 
            plot = FALSE)$amp
        specA_right <- spectro(right, f = samplingrate, wl = wlen, 
            plot = FALSE)$amp
        rm(left, right)
        if (max_freq > nyquist_freq) {
            cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", 
                nyquist_freq, "Hz. But the script was set to measure up to ", 
                max_freq, "Hz.\n\n", sep = ""))
            break
        }
        Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)
        Score <- rep(NA, length(Freq))
        for (j in 1:length(Freq)) {
            Score[j] = getscore(specA_left, Freq[j], (Freq[j] + 
                freq_step), db_threshold, freq_per_row)
        }
        left_vals = Score
        Score1 = 0
        for (i in 1:length(Freq)) {
            Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
        }
        Score_left = (-(Score1))/length(Freq)
        Score <- rep(NA, length(Freq))
        for (j in 1:length(Freq)) {
            Score[j] = getscore(specA_right, Freq[j], (Freq[j] + 
                freq_step), db_threshold, freq_per_row)
        }
        right_vals = Score
        Score1 = 0
        for (i in 1:length(Freq)) {
            Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
        }
        Score_right = (-(Score1))/length(Freq)
        cat(" ==============================================\n")
        cat(paste(" Results for file '", filename, "' \n  with a dB threshold of ", 
            db_threshold, "\n\n", sep = ""))
        cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
        cat("Frequency range (Hz), left channel proportion, right channel proportion\n")
        for (j in seq(length(Freq), 1, by = -1)) {
            cat(paste(Freq[j], "-", (Freq[j] + freq_step), ",", 
                round(left_vals[j], 6), ",", round(right_vals[j], 
                  6), "\n", sep = ""))
        }
        cat("\n Plot of proportions in each band: \n\n")
        cat("  Left channel\n")
        cat("   Freq. range (Hz) \t |--------------------|\n")
        for (j in seq(length(Freq), 1, by = -1)) {
            cat(paste("   ", Freq[j], "-", (Freq[j] + freq_step), 
                " \t\t |", sep = ""))
            temp_val = round(left_vals[j], 2) * 20
            if (temp_val > 0) {
                for (i in 1:temp_val) {
                  cat("*")
                }
            }
            cat("\n")
            rm(temp_val)
        }
        cat("\n  Right channel\n")
        cat("   Freq. range (Hz) \t |--------------------|\n")
        for (j in seq(length(Freq), 1, by = -1)) {
            cat(paste("   ", Freq[j], "-", (Freq[j] + freq_step), 
                " \t\t |", sep = ""))
            temp_val = round(right_vals[j], 2) * 20
            if (temp_val > 0) {
                for (i in 1:temp_val) {
                  cat("*")
                }
            }
            cat("\n")
            rm(temp_val)
        }
        cat("\n  Acoustic diversity (Shannon's Index): \n")
        cat(paste("   Left channel: ", round(Score_left, 6), 
            "\n", sep = ""))
        cat(paste("   Right channel: ", round(Score_right, 6), 
            "\n\n", sep = ""))
        left_adi_return = round(Score_left, 6)
        right_adi_return = round(Score_right, 6)
        cat("  Band Eveness (Gini coefficient): \n")
        cat(paste("   Left channel: ", round(Gini(left_vals), 
            6), "\n", sep = ""))
        cat(paste("   Right channel: ", round(Gini(right_vals), 
            6), "\n\n", sep = ""))
        left_gini_return = round(Gini(left_vals), 6)
        right_gini_return = round(Gini(right_vals), 6)
    }
    else {
        cat("\n This is a mono file.\n")
        cat("\n Getting values from spectrogram... Please wait... \n")
        specA_left <- spectro(soundfile, f = samplingrate, wl = wlen, 
            plot = FALSE)$amp
        rm(soundfile)
        if (max_freq > nyquist_freq) {
            cat(paste("\n ERROR: The maximum acoustic frequency that this file can use is ", 
                nyquist_freq, "Hz. But the script was set to measure up to ", 
                max_freq, "Hz.\n\n", sep = ""))
            break
        }
        Freq <- seq(from = 0, to = max_freq - freq_step, by = freq_step)
        Score <- rep(NA, length(Freq))
        for (j in 1:length(Freq)) {
            Score[j] = getscore(specA_left, Freq[j], (Freq[j] + 
                freq_step), db_threshold, freq_per_row)
        }
        left_vals = Score
        Score1 = 0
        for (i in 1:length(Freq)) {
            Score1 = Score1 + (Score[i] * log(Score[i] + 1e-07))
        }
        Score_left = (-(Score1))/length(Freq)
        cat(" ==============================================\n")
        cat(paste(" Results for file '", filename, "' \n  with a dB threshold of ", 
            db_threshold, "\n\n", sep = ""))
        cat(" Proportion over threshold for each frequency band (in csv format): \n\n")
        cat("Frequency range (Hz), proportion\n")
        for (j in seq(length(Freq), 1, by = -1)) {
            cat(paste(Freq[j], "-", (Freq[j] + freq_step), ",", 
                round(left_vals[j], 6), "\n", sep = ""))
        }
        cat("\n Plot of proportions in each band: \n\n")
        cat("   Freq. range (Hz) \t |--------------------|\n")
        for (j in seq(length(Freq), 1, by = -1)) {
            cat(paste("   ", Freq[j], "-", (Freq[j] + freq_step), 
                " \t\t |", sep = ""))
            temp_val = round(left_vals[j], 2) * 20
            if (temp_val > 0) {
                for (i in 1:temp_val) {
                  cat("*")
                }
            }
            cat("\n")
            rm(temp_val)
        }
        cat("\n  Acoustic diversity (Shannon's Index): ")
        cat(paste(round(Score_left, 6), "\n", sep = ""))
        left_adi_return = round(Score_left, 6)
        right_adi_return = 0
        cat("  Band Eveness (Gini coefficient): ")
        cat(paste(round(Gini(left_vals), 6), "\n", sep = ""))
        left_gini_return = round(Gini(left_vals), 6)
        right_gini_return = 0
    }
    return(c(left_adi_return, right_adi_return, left_gini_return, 
        right_gini_return))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
