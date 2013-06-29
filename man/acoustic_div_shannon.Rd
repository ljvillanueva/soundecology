\name{acoustic_diversity}

\alias{acoustic_diversity}

\title{Acoustic Diversity Index}

\description{Acoustic Complexity Index from Pieretti, \emph{et al.} 2011.}

\usage{
acoustic_diversity(soundfile, save_to_text=FALSE, max_freq=10000, db_threshold=-50, freq_step=1000)
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
 Pieretti, N., A. Farina, and D. Morri. 2011. A new methodology to infer the singing activity of an avian community: The Acoustic Complexity Index (ACI). Ecological Indicators 11: 868-873.}

\keyword{soundecol}
\keyword{acoustic_diversity}