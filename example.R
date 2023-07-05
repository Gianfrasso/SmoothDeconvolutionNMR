# Libraries
rm(list = ls())
source("support_functions_nmr.R")

# Load the data
data = readMat("NMRpotato/go01.mat")
ix = 33
Fx = (data$X[ix, ])
x = c(data$Time)
n = length(x)

# Smooth deconvolution
ests = smooth_deconvolution(Fx = Fx)
h = ests$h
la = ests$la
resid = ests$resid
mu = ests$mu 
res = ests$res
lla = ests$lla
tt = ests$tt 
xa = ests$xa

# Plot results
fitDat = data.frame(
  x = x, obs = Fx, fit = mu,
  resid = res 
)
estDat = data.frame(time = 10^tt, h = h)
penDat = data.frame(
  xa = xa, pen = lla,
  ta = 10^seq(-3.5, 3.5, len = length(xa))
)

fitPlt = ggplot(fitDat, aes(x = x, y = obs)) +
  geom_point() +
  geom_line(aes(x = x, y = fit), colour = "orange") +
  xlab("Time (ms)") +
  ylab("Residual magnetisation") +
  ggtitle(paste("Sample", ix, "(DM =", paste0(round(data$Y[ix, 1], 5), "%)"))) +
  theme_bw()

estPlt = ggplot(estDat, aes(x = time, y = h)) +
  geom_line() +
  scale_x_continuous(trans = "log10", 
                    limits = range(data$Time) + c(-0.6, 1.3e3)) + 
  xlab("Time (ms) - log10 scale") +
  ylab("Density") +
  theme_bw()

resPlt = ggplot(fitDat, aes(x = x, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0.0, colour = "gray50") +
  xlab("Time (ms)") +
  ylab("Residuals") +
  theme_bw()

penPlt = ggplot(penDat, aes(x = ta, y = pen)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(trans = "log10") +
  xlab("Time (ms) - log10 scale") +
  ylab("Adaptive penalty - log10 scale") +
  theme_bw()

plRes = gridExtra::grid.arrange(fitPlt, estPlt,
  resPlt, penPlt,
  nrow = 2, ncol = 2
)
ggsave("Figure2.pdf",plRes, width = 10, height = 6) #Example
