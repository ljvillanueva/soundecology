#Soundscape Index test
#
#
#Gage's Soundscape Index
# From REAL - Remote Environmental Assessment Laboratory
# http://www.real.msu.edu/


require("tuneR")
require("seewave")
require("pracma")
#require("ineq")

source("R/soundscape_index.R")

soundfile <- readWave("test_scripts/testsound.wav")

vals <- ndsi(soundfile)
