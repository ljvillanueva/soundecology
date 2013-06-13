#Acoustic diversity index test
#

require("tuneR")
require("seewave")
#require("pracma")
#require("ineq")

source("R/seewave_simplified.R")
source("R/acoustic_complexity_index.R")

soundfile <- readWave("test_scripts/testsound.wav")

vals <- acoustic_complexity(soundfile)
