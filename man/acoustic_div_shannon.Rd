\name{acoustic_diversity}

\alias{acoustic_diversity}

\title{Acoustic Diversity Index}

\description{Acoustic Diversity Index from Villanueva-Rivera \emph{et al.} 2011.}

\usage{
acoustic_diversity(soundfile, save_to_text=FALSE, max_freq=10000, db_threshold=-50, freq_step=1000)
}

\arguments{
	\item{soundfile}{ an object of class \code{"Wave"} loaded with the function readWave of the tuneR package. }
	\item{save_to_text}{ boolean to save the values to a text file. }
	\item{max_freq}{ maximum frequency to use when calculating the value, in Hertz. }
	\item{db_threshold}{ threshold to use in the calculation. }
	\item{freq_step}{ size of frequency bands. }
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
Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.}

\keyword{soundecol}
\keyword{acoustic_diversity}