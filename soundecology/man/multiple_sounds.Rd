\name{multiple_sounds}

\alias{multiple_sounds}

\title{Multiple sound files}

\description{Function to extract the specified index from all the wav files in a directory.}

\usage{
multiple_sounds(directory, resultfile, soundindex, ...)
}

\arguments{
	\item{directory}{ a valid directory, readable by the user, that contains the wav files. }
	\item{resultfile}{ name of the text file to which write the results in csv format.}
	\item{soundindex}{ which index to calculate: 
		\itemize{
         \item \code{ndsi}
         \item \code{acoustic_complexity}
         \item \code{acoustic_diversity}
         \item \code{bioacoustic_index}
         \item \code{H} from the \code{seewave} package
       }
	}
	\item{...}{ additional variables to pass to the selected function. See each function's help for details. }
	}

\author{Luis J. Villanueva-Rivera}

\keyword{soundecology}
\keyword{multiple_sounds}
