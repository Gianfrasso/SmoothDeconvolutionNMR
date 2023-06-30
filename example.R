# Libraries
rm(list = ls())
source("support_functions_nmr.R")

# Load the data
data = readMat("NMRpotato/go01.mat")
ix = 33
scl = max(data$X[ix, ])
Fx = (data$X[ix, ]) / scl
x = c(data$Time)
n = length(x)

# Define model matrix (with integration rule)
min_t = -3.5
max_t = 3.5
tt = define_t2_vector(m = 200, min_t = min_t, max_t = max_t)
Ci = define_model_matrix(x, tt)

# Define bases and penalty
bdeg = 3
ndx = 35
dd = 2
B = as.spam(bbase(tt, bdeg = bdeg, nseg = ndx))
nb = ncol(B)
D = diff(diag(nb), diff = dd)

# Bases for adaptive
xa = 1:(nb - dd) / (nb - dd)
Ba = bbase(xa, nseg = 5, bdeg = 3)
pen_comp = define_adaptive_penalty(Ba, D, dd, extra_ridge=1e-14)
Pl = pen_comp$Pl
Ll = pen_comp$Ll

# Estimate eta/lambdas
ests = fit_smooth_deconvolution(Fx = Fx, Ci = Ci, Pl = Pl, B = B)
h = ests$h * scl
la = ests$la
resid = ests$resid
mu = Ci %*% h 
res = ests$res * scl

# Plot results
fitDat = data.frame(
  x = x, obs = Fx * scl, fit = mu,
  resid = res 
)
estDat = data.frame(time = 10^tt, h = h)
penDat = data.frame(
  xa = xa, pen = log10(Ba %*% la),
  ta = 10^seq(min_t, max_t, len = nrow(Ba))
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
