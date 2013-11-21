\name{acoustic_evenness}

\alias{acoustic_evenness}

\title{Acoustic Evenness Index}

\description{Acoustic Evenness Index from Villanueva-Rivera \emph{et al.} 2011 (band evenness using the Gini index). 
The AEI is calculated by dividing the spectrogram into bins (default 10) and taking the proportion of the signals in each bin above a threshold (default -50 dBFS). The AEI is the result of the Gini index applied to these bins.}

\usage{
acoustic_evenness(soundfile, max_freq = 10000, db_threshold = -50, freq_step = 1000)
}

\arguments{
	\item{soundfile}{ an object of class \code{Wave} loaded with the function readWave of the \code{tuneR} package. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{db_threshold}{ threshold to use in the calculation. }
	\item{freq_step}{ size of frequency bands. }
	}

\value{
	Returns a list with five objects per channel
	\item{aei_left}{ AEI for the left channel }
	\item{aei_right}{ AEI for the right channel }
	}

\note{For audio files with one channel, the results are showed as the left channel, the right channel returns \code{NA}.} 
  
\examples{
	data(tropicalsound)
	result <- acoustic_evenness(tropicalsound)
	print(result$aei_left)
	
	summary(result)
	}

\references{
Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.}

\keyword{soundecology}
\keyword{acoustic_evenness}
