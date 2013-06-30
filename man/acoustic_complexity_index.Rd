\name{acoustic_complexity}

\alias{acoustic_complexity}

\title{Acoustic Complexity Index}

\description{Acoustic Complexity Index from Pieretti, \emph{et al.} 2011. 
The ACI is based on the "observation that many biotic sounds, such as bird songs, are 
characterized by an intrinsic variability of intensities, while some types
of human generated noise (such as car passing or airplane transit)
present very constant intensity values" (Pieretti, \emph{et al.} 2011).

The index was found to be equivalent to the one calculated using SoundscapeMeter v1.0.14.05.2012, courtesy of A. Farina.}

\usage{
	acoustic_complexity(soundfile, max_freq=22050, j=5, fft_w=512)
	}

\arguments{
	\item{soundfile}{ an object of class Wave loaded with the function readWave of the tuneR package. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{j}{ the cluster size, in seconds. }
	\item{fft_w}{ FFT window to use. }
	}

\value{
	Returns a list with three objects per channel
	\item{AciTotAll_left}{ the ACI total for the left channel }
	\item{AciTotAll_right}{ the ACI total for the right channel }
	\item{aci_fl_left_vals}{ values of ACI(fl) for the left channel }
	\item{aci_fl_right_vals}{ values of ACI(fl) for the right channel }
	\item{aci_left_matrix}{ Matrix of all values before calculating ACI(fl) for the left channel }
	\item{aci_right_matrix}{ Matrix of all values before calculating ACI(fl) for the right channel }
	}

\author{Luis J. Villanueva-Rivera}

\examples{
	data(tropicalsound)
	ACI <- acoustic_complexity(tropicalsound)
	summary(ACI)
	print(ACI$AciTotAll_left)
	}

\references{
 Pieretti, N., A. Farina, and D. Morri. 2011. A new methodology to infer the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators 11: 868-873. doi: 10.1016/j.ecolind.2010.11.005}

\keyword{soundecology}
\keyword{acoustic_complexity}
