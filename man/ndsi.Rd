\name{ndsi}

\alias{ndsi}

\title{Normalized Difference Soundscape Index}

\description{Normalized Difference Soundscape Index from REAL and Krause, \emph{et al.} 2011.}

\usage{
ndsi(soundfile, anthro_min=1000, anthro_max=2000, bio_min=2000, bio_max=8000, hz_interval=1000)
}

\arguments{
	\item{soundfile}{ an object of class \code{"Wave"} loaded with the function readWave of the tuneR package. }
	\item{anthro_min}{ minimum value of the range of frequencies of the anthrophony.}
	\item{anthro_max}{ maximum value of the range of frequencies of the anthrophony.}
	\item{bio_min}{ minimum value of the range of frequencies of the biophony.}
	\item{bio_max}{ maximum value of the range of frequencies of the biophony.}
	\item{hz_interval}{ frequency intervals used in the calculation, in Hertz. }
	}

\value{
	Returns a list with three objects per channel
	\item{ndsi_left}{ NDSI value for the left channel }
	\item{ndsi_right}{ NDSI value for the left channel }
	}

\author{Luis J. Villanueva-Rivera}

%\examples{
	%library("seewave")
	%data(tropicalsound)
	%ACI <- acoustic_complexity(tropicalsound)
	%summary(ACI)
	%print(ACI$AciTotAll_left)
	%}

\references{
Remote Environmental Assessment Laboratory. http://www.real.msu.edu

Krause, Bernie, Stuart H. Gage, and Wooyeong Joo. 2011. Measuring and interpreting the temporal variability in the soundscape at four places in Sequoia National Park. Landscape Ecology 26: 1247-1256. doi: 10.1007/s10980-011-9639-6.}

\keyword{soundecol}
\keyword{ndsi}