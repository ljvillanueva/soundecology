\name{multiple_sounds}

\alias{multiple_sounds}

\title{Multiple sound files}

\description{Function to extract the specified index from all the wav or flac files in a directory. The results, including the filename and wave technical details, are saved to a csv file. If the computer has multiple cores, it can run files in parallel.}

\usage{
multiple_sounds(directory, resultfile, soundindex, no_cores=1, flac = FALSE, ...)
}

\arguments{
	\item{directory}{ a valid directory, readable by the user, that contains the wav files. }
	\item{resultfile}{ name of the text file to which write the results in comma-separated values format.}
	\item{soundindex}{ which index to calculate: 
		\itemize{
         \item \code{ndsi}
         \item \code{acoustic_complexity}
         \item \code{acoustic_diversity}
         \item \code{acoustic_eveness}
         \item \code{bioacoustic_index}
         \item \code{H} from the \code{seewave} package
       }
	}
	\item{no_cores}{ number of cores to use when calculating the indices. Can be \code{max} to use all cores, \code{-1} to use all but one core, or any positive integer. Default is 1. Uses the \code{parallel} package.}
	\item{flac}{ logical variable to indicate that the files are in flac format. Flac must be installed in the system. Uses the function \code{wav2flac} of \code{seewave}.}
	\item{...}{ additional variables to pass to the selected function. See each function's help for details.}
	}

\examples{
	\dontrun{
		#Calculate the ACI of all the wav
		# files in the directory "/home/user/wavs/"
		# using the function acoustic_complexity
		multiple_sounds(directory="/home/user/wavs/", 
			resultfile="/home/user/results.csv", 
			soundindex="acoustic_complexity")
		
		#Calculate the same as above using 12000Hz as the
		# maximum frequency instead of the default.
		multiple_sounds(directory="/home/user/wavs/", 
			resultfile="/home/user/results.csv", 
			soundindex="acoustic_complexity", max_freq=12000)
			
		#Calculate the same as above using two cores
		multiple_sounds(directory="/home/user/wavs/", 
			resultfile="/home/user/results.csv", 
			soundindex="acoustic_complexity", no_cores=2)
			
		#Calculate the same as above using all the cores
		# the computer has
		multiple_sounds(directory="/home/user/wavs/", 
			resultfile="/home/user/results.csv", 
			soundindex="acoustic_complexity", no_cores="max")
			
		#Calculate the same as above using all but one cores
		multiple_sounds(directory="/home/user/wavs/", 
			resultfile="/home/user/results.csv", 
			soundindex="acoustic_complexity", no_cores=-1)
		}
	}

\keyword{soundecology}
\keyword{multiple_sounds}
