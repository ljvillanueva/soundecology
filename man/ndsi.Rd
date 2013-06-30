\name{ndsi}

\alias{ndsi}

\title{Normalized Difference Soundscape Index}

\description{Normalized Difference Soundscape Index from REAL and Kasten, \emph{et al.} 2012.}

\usage{
ndsi(soundfile, anthro_min=1000, anthro_max=2000, 
	bio_min=2000, bio_max=8000, hz_interval=1000)
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

Kasten, Eric P., Stuart H. Gage, Jordan Fox, and Wooyeong Joo. 2012.
The Remote Environmental Assessment Laboratory’s Acoustic Library: An Archive for 
Studying Soundscape Ecology. Ecological Informatics 12: 50-67. doi: 10.1016/j.ecoinf.2012.08.001}

\keyword{soundecology}
\keyword{ndsi}
