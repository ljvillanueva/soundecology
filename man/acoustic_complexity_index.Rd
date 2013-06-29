\name{acoustic_complexity}

\alias{acoustic_complexity}

\title{Acoustic Complexity Index}

\description{Acoustic Complexity Index from Pieretti, \emph{et al.} 2011.}

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
	An object of class \code{Wave}, basically a list including elements
	\item{coefficients}{ a named vector of coefficients }
	\item{vcov}{ covariance matrix of coefficients }
	\item{fitted.values}{ fitted values }
	\item{residuals}{ residuals }
	}

\author{Luis J. Villanueva-Rivera}

\examples{
	library("seewave")
	data(tropicalsound)
	ACI <- acoustic_complexity(tropicalsound)
	summary(ACI)
	print(ACI$AciTotAll_left)
	}

\references{
 Pieretti, N., A. Farina, and D. Morri. 2011. A new methodology to infer the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators 11: 868-873. doi: 10.1016/j.ecolind.2010.11.005}

\keyword{soundecol}
\keyword{acoustic_complexity}