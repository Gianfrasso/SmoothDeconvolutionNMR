# Libraries
rm(list = ls())
source("smooth_deconvolution.R")

# Load the data
data = readMat("NMRpotato/go04.mat")
Fx = c(data$X)
x = c(data$Time)
n = length(x)

# Smooth deconvolution
ests = smooth_deconvolution(Fx = Fx, x = x)
plRes = plot(ests)
ggsave("Figure_example2.pdf", plRes, width = 10, height = 6)  #Example
