\name{bioacoustic_index}

\alias{bioacoustic_index}

\title{Bioacoustic Index}

\description{Bioacoustic Index from Boelman, \emph{et al.} 2007. 
Inspired on Matlab code courtesy of NT Boelman.
Several parts where changed, in particular log math, so this won't be
directly comparable to the original code in the paper.

The Bioacoustic Index is calculated as the "area under each curve included all frequency bands associated
with the dB value that was greater than the minimum dB value for each curve. The area values are thus a
function of both the sound level and the number of frequency bands used by the avifauna" (Boelman, \emph{et al.} 2007).}

\usage{
bioacoustic_index(soundfile, min_freq = 2000, max_freq = 8000, fft_w = 512)
}

\arguments{
	\item{soundfile}{ an object of class \code{Wave} loaded with the function readWave of the \code{tuneR} package. }
	\item{min_freq}{ minimum frequency to use when calculating the value, in Hertz. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{fft_w}{ FFT window size. }
	}

\value{
	Returns a list with one object per channel
	\item{left_area}{ area under the curve for the left channel }
	\item{right_area}{ area under the curve for the right channel }
	}

\examples{
	data(tropicalsound)
	bioindex <- bioacoustic_index(tropicalsound)
	print(bioindex$left_area)
	
	summary(bioindex)
	}

\references{
Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecological Applications 17: 2137-2144.}

\keyword{soundecology}
\keyword{acoustic_diversity}
