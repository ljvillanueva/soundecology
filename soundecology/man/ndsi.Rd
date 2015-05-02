\name{ndsi}

\alias{ndsi}

\title{Normalized Difference Soundscape Index}

\description{Normalized Difference Soundscape Index (NDSI) from REAL and Kasten, \emph{et al.} 2012. The NDSI
seeks to "estimate the level of anthropogenic disturbance on the soundscape by
computing the ratio of human-generated (anthrophony) to biological
(biophony) acoustic components found in field collected sound samples" (Kasten, \emph{et al.} 2012).

Tested with Matlab code courtesy of S. Gage.}

\usage{
ndsi(soundfile, fft_w = 1024, anthro_min = 1000, anthro_max = 2000, 
	bio_min = 2000, bio_max = 11000)
}

\arguments{
	\item{soundfile}{ an object of class \code{Wave} loaded with the function readWave of the \code{tuneR} package. }
	\item{fft_w}{ FFT window size.}
	\item{anthro_min}{ minimum value of the range of frequencies of the anthrophony.}
	\item{anthro_max}{ maximum value of the range of frequencies of the anthrophony.}
	\item{bio_min}{ minimum value of the range of frequencies of the biophony.}
	\item{bio_max}{ maximum value of the range of frequencies of the biophony.}
	}

\value{
	Returns a list with one object per channel
  \item{ndsi_left}{ NDSI value for the left channel }
	\item{ndsi_right}{ NDSI value for the right channel }
  \item{biophony_left}{ value for the biophony for the left channel }
	\item{anthrophony_left}{ value for the anthrophony for the left channel }
  \item{biophony_right}{ value for the biophony for the right channel }
	\item{anthrophony_right}{ value for the anthrophony for the right channel }
	}

\details{The bin size is determined as the difference between anthro_max and anthro_min, by default 1000 Hz.
}  


\examples{
	data(tropicalsound)
	NDSI <- ndsi(tropicalsound)
	print(NDSI$ndsi_left)
	
	summary(NDSI)
	}

\references{
Remote Environmental Assessment Laboratory. http://www.real.msu.edu

Kasten, Eric P., Stuart H. Gage, Jordan Fox, and Wooyeong Joo. 2012.
The Remote Environmental Assessment Laboratory's Acoustic Library: An Archive for 
Studying Soundscape Ecology. Ecological Informatics 12: 50-67. doi: 10.1016/j.ecoinf.2012.08.001
}

\keyword{soundecology}
\keyword{ndsi}
