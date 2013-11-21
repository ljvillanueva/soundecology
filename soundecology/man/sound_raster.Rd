\name{sound_raster}

\alias{sound_raster}

\title{ASCII raster from sound file}

\description{This function creates a raster file in ASCII format from the spectrogram of a soundfile. This file can be opened in ArcGIS or any other GIS software. For more details see the tutorial of Villanueva-Rivera \emph{et al.} 2011.}

\usage{
sound_raster(wavfile = NA, wav_directory = NA, max_freq = 10000, no_cores = 1)
}

\arguments{
	\item{wavfile}{ a single sound file in wav format. }
	\item{max_freq}{ maximum frequency to draw the spectrogram, in Hertz. }
	\item{wav_directory}{ a directory that contains wav files. To specify the working directory, use \code{wav_directory="."} }
	\item{no_cores}{ number of cores to use when working in a directory. Can be \code{max} to use all cores, \code{-1} to use all but one core, or any positive integer. Default is 1. Uses the \code{parallel} package.}
	}

\value{The function will save a file for each channel, in the same directory where the files are at, with the extension .asc.}

\note{ To get a raster file for a single file, use the argument \code{wavfile}. For many files, use the argument \code{wav_directory}. Do not use both at the same time or the function will return an error. 

This function was released with the version 1.3 of the tutorial of the primer paper, available at:

http://ltm.agriculture.purdue.edu/soundscapes/primer/

and at the website of the package:

http://ljvillanueva.github.io/soundecology/}

\examples{
\dontrun{
sound_raster(wavfile = "file1.wav")

sound_raster(wav_directory = "/home/user/wavdirectory")

sound_raster(wav_directory = "/home/user/wavdirectory", no_cores = 4)
	}
	}

\references{
Villanueva-Rivera, L. J., B. C. Pijanowski, J. Doucette, and B. Pekin. 2011. A primer of acoustic analysis for landscape ecologists. Landscape Ecology 26: 1233-1246. doi: 10.1007/s10980-011-9636-9.}

\keyword{soundecology}
\keyword{sound_raster}
