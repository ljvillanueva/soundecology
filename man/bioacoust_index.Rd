\name{bioacoustic_index}

\alias{bioacoustic_index}

\title{Bioacoustic Index}

\description{Bioacoustic Index from Boelman, \emph{et al.} 2007. 
Inspired on Matlab code courtesy of NT Boelman.}

\usage{
bioacoustic_index(soundfile, min_freq=2000, max_freq=8000, fft_w=512)
}

\arguments{
	\item{soundfile}{ an object of class \code{"Wave"} loaded with the function readWave of the tuneR package. }
	\item{min_freq}{ minimum frequency to use when calculating the value, in Hertz. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{fft_w}{ FFT window size. }
	}

\value{
	Returns a list with three objects per channel
	\item{left_area}{ area under the curve for the left channel }
	\item{right_area}{ area under the curve for the right channel }
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
Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17:2137-2144.}

\keyword{soundecology}
\keyword{acoustic_diversity}
