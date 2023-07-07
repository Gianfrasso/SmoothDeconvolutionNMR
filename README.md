# SmoothDeconvolutionNMR
This repository contains the main functionalities to resolve peaks in NMR relaxation experiments using the smooth deconvolution approach of Frasso and Eilers (2023, forthcoming).
The main computations are performed by `smooth_deconvolution()` function included in the `support_functions_nmr.R` module.

The repository also contains the potato tubers data presented in Hansen, C. L., Thybo, A. K., Bertram, H. C., Viereck, N., van den Berg, F., and Engelsen, S. B. (2010). 
Determination of dry matter content in potato tubers by low-field nuclear magnetic resonance (LF-NMR). _Journal of Agricultural and Food Chemistry_, 58(19), 10300–10304.

Running the example below (see also `example.R`) file it is possible to reproduce Figure 2 of the reference.
A package is under construction. In the meantime, please contact me at gianfrasso@gmail.com for questions and issues.

```r
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
ggsave("Figure2.pdf", plRes, width = 10, height = 6)  
```
