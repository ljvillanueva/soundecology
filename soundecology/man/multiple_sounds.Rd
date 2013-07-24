\name{multiple_sounds}

\alias{multiple_sounds}

\title{Multiple sound files}

\description{Function to extract the specified index from all the wav files in a directory.}

\usage{
multiple_sounds(directory, filename, soundindex, ...)
}

\arguments{
	\item{directory}{ a valid directory, readable by the user, that contains the wav files. }
	\item{filename}{ name of the text file to which write the results in csv format.}
	\item{soundindex}{ which index to calculate: \code{"ndsi"}, \code{"acoustic_complexity"}, \code{"acoustic_diversity"}, or \code{"bioacoustic_index"}.}
	\item{...}{ additional variables to give the appropriate function. See each function's help for details. }
	}

\author{Luis J. Villanueva-Rivera}

\keyword{soundecology}
\keyword{multiple_sounds}
