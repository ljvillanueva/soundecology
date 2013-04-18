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

source("R/seewave_simplified.R")

soundfile1 <- readWave("test_scripts/testsound.wav")
soundfile2 <- readWave("test_scripts/testsound2.wav")

spec1 <- meanspec(soundfile1, plot=FALSE)
spec2 <- meanspec(soundfile2, plot=FALSE)

correlation <- corspec(spec1, spec2, plot=FALSE)
        
correlation$rmax
correlation$p
correlation$f



#diffwave
if (soundfile1@samp.rate == soundfile2@samp.rate){
  
  min_length <- min(length(soundfile1@left)/soundfile1@samp.rate, length(soundfile2@left)/soundfile1@samp.rate)
  
  soundfile1_1 <- cutw(soundfile1, from=0, to=min_length)
  soundfile2_1 <- cutw(soundfile2, from=0, to=min_length)
  
  this_diff <- diffwave(soundfile1_1, soundfile2_1, f=soundfile1@samp.rate)
  cat("\n diffwave: ")
  cat(this_diff)
}else{
  cat("\n Waves do not match in sampling rate.")
} 


#Q
Q(meanspec(soundfile1, dB="max0"))


#roughness(spec)
roughness(meanspec(soundfile1, plot=FALSE))



#rugo
spec_1 <- meanspec(soundfile1, plot=FALSE)[,2]
rugo(spec_1/max(spec_1))



#sfm
spec_1 <- meanspec(soundfile1, plot=FALSE)[,2]
sfm(spec(soundfile2, plot=FALSE)[,2])



#sh
sh(spec(soundfile2, plot=FALSE)[,2])



#simspec
simspec(spec(soundfile1_1, f=soundfile1@samp.rate), spec(soundfile2_1, f=soundfile1@samp.rate), f=soundfile1@samp.rate)
