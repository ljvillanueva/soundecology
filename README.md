[![Build Status](https://travis-ci.org/ljvillanueva/soundecology.svg?branch=master)](https://travis-ci.org/ljvillanueva/soundecology)

soundecology
=========

R package with functions to calculate indices for soundscape ecology and other ecology research that uses audio recordings.

Indices:

 * Normalized Difference Soundscape Index (NDSI) from REAL and Kasten et al. 2012.
 * Acoustic Complexity Index (ACI) from Pieretti et al. 2011. 
 * Acoustic Diversity Index (ADI) from Villanueva-Rivera et al. 2011.
 * Bioacoustic Index from Boelman et al. 2007. 

Other functions:

 * measure_signals() - This function lets the user select bounding boxes to get statistics of the signals of interest in a sound
file.
 * multiple_sounds() - Function to extract the specified index from all the wav or flac files in a directory. The results,
including the filename and wave technical details, are saved to a csv file. If the computer has
multiple cores, it can run files in parallel.
 * sound_raster() - This function creates a raster file in ASCII format from the spectrogram of a soundfile.

The package is available from CRAN: http://cran.r-project.org/web/packages/soundecology/

http://ljvillanueva.github.io/soundecology/
