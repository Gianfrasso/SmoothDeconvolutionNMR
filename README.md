# SmoothDeconvolutionNMR
This repository contains the main functionalities to resolve peaks in NMR relaxation experiments using the smooth deconvolution approach of Gianluca Frasso, Paul H.C. Eilers (2024),
Smooth deconvolution of low-field NMR signals,
_Analytica Chimica Acta_ (https://doi.org/10.1016/j.aca.2023.341808).
The main computations are performed by `smooth_deconvolution()` function included in the `methods.R` module.

The repository also contains the potato tubers data presented in Hansen, C. L., Thybo, A. K., Bertram, H. C., Viereck, N., van den Berg, F., and Engelsen, S. B. (2010). 
Determination of dry matter content in potato tubers by low-field nuclear magnetic resonance (LF-NMR). _Journal of Agricultural and Food Chemistry_, 58(19), 10300–10304.

## Installation

1. Clone the repository.
2. Open R in the project directory.
3. Install `renv` (if you don't have it already): `install.packages("renv")`
4. Restore the project's library: `renv::restore()`

This will install the exact versions of the R packages used in this project, as specified in the `renv.lock` file.

## Usage

After restoring the environment, you can run the R scripts in the project. 

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

# Smooth deconvolution
ests = smooth_deconvolution(Fx = Fx, x = x)
plRes = plot(ests)
ggsave("Figure2.pdf", plRes, width = 10, height = 6)  
```

