\name{measure_signals}

\alias{measure_signals}

\title{Measure a signal or song in a wavefile}

\description{This function lets the user select bounding boxes to get statistics of the signals of interest in a sound file. }

\usage{
measure_signals(wavfile, wl = 512, min_freq = NA, max_freq = NA, min_time = NA,
  max_time = NA, plot_range = 50, dBFS_range = 30, sample_size = 1,
  resultfile = NA, channel = "left")
}

\arguments{
	\item{wavfile}{ a sound file in wav format. }
	\item{wl}{ window length for the spectrogram. }
	\item{min_freq}{ minimum frequency to draw the spectrogram, in kiloHertz. }
	\item{max_freq}{ maximum frequency to draw the spectrogram, in kiloHertz. }
	\item{min_time}{ minimum time to draw the spectrogram, in seconds. }
	\item{max_time}{ maximum time to draw the spectrogram, in seconds. }
	\item{plot_range}{ lower limit of values to plot the spectrogram. }
	\item{dBFS_range}{ range of values that is considered a signal, based on the maximum that is calculated. See notes below. }
	\item{sample_size}{ number of samples to measure in the spectrogram. }
	\item{resultfile}{ name of the file to save the results. }
	\item{channel}{ which channel to plot. }
	}

\value{The function will open a spectrogram plot to allow the user to click on the regions of interest. Once all the samples are selected, the function saves a file with the values measured in each sample. In addition, the results of the function \code{dfreq} of the package \code{seewave} are saved on a folder named the same as the wavfile, without the .wav extension. }

\note{ For the \code{dBFS_range} argument, the code uses the maximum of the values inside the selected region and saves as a resulting signal the values that fall between \code{(maximum - dBFS_range)} and the maximum. A selected region with a maximum value of -5 and \code{dBFS_range} set to 30 will consider the area with values between -35 and -5 dBFS as a signal. 

The function creates a folder \code{dfreq} where it saves csv files with the results of the function \code{dfreq} from \code{seewave}. The name of each file is coded as: wavfile.samplenumber.csv
}

\examples{
\dontrun{
#Take 5 samples of the file file.wav between 1 - 4 kHz, from 10 to 30 seconds.
measure_signals(wavfile="file.wav", wl=2048, min_freq=1, max_freq=4,
  dBFS_range=30, min_time=10, max_time=30, sample_size=5,
  resultfile="results.csv", plot_range=70)
	}
	}

\keyword{soundecology}
\keyword{measure_signals}