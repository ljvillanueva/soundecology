\name{bioacoustic_index}

\alias{bioacoustic_index}

\title{Bioacoustic Index}

\description{Bioacoustic Index from Boelman, \emph{et al.} 2007. 
Inspired on Matlab code courtesy of NT Boelman.}

\usage{
bioacoustic_index(soundfile, save_to_text=FALSE, max_freq=10000, db_threshold=-50, freq_step=1000)
}

\arguments{
	\item{soundfile}{ an object of class \code{"Wave"} loaded with the function readWave of the tuneR package. }
	\item{save_to_text}{ boolean to save the values to a text file. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{db_threshold}{ threshold to use in the calculation. }
	\item{freq_step}{ size of frequency bands. }
	\item{\dots}{ not used. }
	}

\value{
	Returns a list of values with the elements:
	\item{coefficients}{ a named vector of coefficients }
	\item{vcov}{ covariance matrix of coefficients }
	\item{fitted.values}{ fitted values }
	\item{residuals}{ residuals }
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
Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion resistance in Hawaii: bioacoustics, field surveys, and airborne remote sensing. Ecol Applications 17(8):2137-44.}

\keyword{soundecol}
\keyword{acoustic_diversity}