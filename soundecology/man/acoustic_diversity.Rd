\name{acoustic_diversity}

\alias{acoustic_diversity}

\title{Acoustic Diversity Index}

\description{Acoustic Diversity Index from Villanueva-Rivera \emph{et al.} 2011. 
The ADI is calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default -50 dBFS). The ADI is the result of the Shannon index applied to these bins.}

\usage{
acoustic_diversity(soundfile, max_freq = 10000, db_threshold = -50, 
freq_step = 1000, shannon = TRUE)
}

\arguments{
	\item{soundfile}{ an object of class \code{Wave} loaded with the function readWave of the \code{tuneR} package. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{db_threshold}{ threshold to use in the calculation. }
	\item{freq_step}{ size of frequency bands. }
	\item{shannon}{ TRUE to use the Shannon's diversity index to calculate the ADI (default). }
	}

\value{
	Returns a list with five objects per channel
	\item{adi_left}{ ADI value for the left channel }
	\item{adi_right}{ ADI value for the right channel }
	\item{left_band_values}{ vector of proportion values for each band for the left channel }
	\item{right_band_values}{ vector of proportion values for each band for the right channel }
	\item{left_bandrange_values}{ vector of frequency values for each band for the left channel }
	\item{right_bandrange_values}{ vector of frequency values for each band for the right channel }
	}

\note{The code to calculate the ADI has changed due to an error we detected in the original scripts in which the value was calculated using a different equation. In a test of ~38k files, both ways to calculate were highly correlated. This version of the function uses the Shannon's Diversity Index. To obtain a result using the old calculation, set the argument \code{shannon} to \code{FALSE}. Please check the vignette "Changes in the Acoustic Diversity Index", included in the package, for more details.

For audio files with one channel, the results are showed as the left channel, the right channel returns \code{NA}.} 

\examples{
	data(tropicalsound)
	result <- acoustic_diversity(tropicalsound)
	%print(result$adis_left)
	print(result$adi_left)
	
	summary(result)
	}

\references{
Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.}

\keyword{soundecology}
\keyword{acoustic_diversity}
