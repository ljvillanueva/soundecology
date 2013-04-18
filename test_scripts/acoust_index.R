#Bioacoustic index
#
#From the paper:
# Boelman NT, Asner GP, Hart PJ, Martin RE. 2007. Multi-trophic invasion
#   resistance in Hawaii: bioacoustics, field surveys, and airborne
#   remote sensing. Ecol Applications 17(8):2137-44.
#
# Inspired on Matlab code provided by NT Boelman. Still need to check if
#  results are comparable.
# Boelman et al. 2007 used min_freq=2000, max_freq=8000, fft_w=1024



require("tuneR")
#require("seewave")
require("pracma")
require("ineq")

source("R/seewave_simplified.R")
source("R/acoustic_div_shannon.R")

soundfile <- readWave("test_scripts/testsound.wav")

acoustic_diversity(soundfile)
