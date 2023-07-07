# Libraries
rm(list = ls())
source("smooth_deconvolution.R")

# Load the data
data = readMat("NMRpotato/go01.mat")
ix = 33
Fx = (data$X[ix, ])
x = c(data$Time)
n = length(x)

# Smooth deconvolution
ests = smooth_deconvolution(Fx = Fx, x = x)
plRes = plot(ests)
ggsave("Figure2.pdf", plRes, width = 10, height = 6)  #Example
