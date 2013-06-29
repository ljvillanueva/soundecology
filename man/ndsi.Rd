\name{ndsi}

\alias{ndsi}

\title{Normalized Difference Soundscape Index}

\description{Normalized Difference Soundscape Index from REAL and Krause, \emph{et al.} 2011.}

\usage{
ndsi(soundfile, save_to_text=FALSE, max_freq=10000, db_threshold=-50, freq_step=1000)
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
Remote Environmental Assessment Laboratory. http://www.real.msu.edu

Krause, Bernie, Stuart H. Gage, and Wooyeong Joo. 2011. Measuring and interpreting the temporal variability in the soundscape at four places in Sequoia National Park. Landscape Ecology 26: 1247-1256. DOI 10.1007/s10980-011-9639-6.}

\keyword{soundecol}
\keyword{ndsi}