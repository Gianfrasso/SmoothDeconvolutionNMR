# Install/load libraries
libs = c("R.matlab", "spam", "JOPS", "scales", "ggplot2")
install.packages(setdiff(libs, rownames(installed.packages())))
pkgLoad = lapply(libs, function(x) require(x, character.only = T))

#' Vector of T2 times
#' 
#' Creates a fector of log10 T2 values
#' @param m: integer length of the T2 times
#' @param min_t: min T2 on log10 scale
#' @param max_t: max T2 on log10 scale
#' @return Vector of T2 values
define_t2_vector = function(m, min_t, max_t) {
    tt = seq(min_t, max_t, len = m)
    return(tt)
}

#' Model matrix
#' 
#' Computes the model matrix for integral approx.
#' @param x: vector of boserved times (dim n)
#' @param tt: vector of m T2 times (log10 scale)
#' @return An integration matrix (n x m)
define_model_matrix = function(x, tt) {
    Ci = exp(outer(-x, 1/(10^tt)))
    return(Ci)
}

#' Adaptive bases penalty
#' 
#' Computes a list of matrices composing the adaptive penalty
#' @param Ba: B-spline basis for adaptive penalty
#' @param D: difference penalty matrix
#' @param dd: integer order of the difference penalty 
#' @param extra_ridge: float to condition P
#' @return list of adaptive penalty components
define_adaptive_penalty = function(Ba, D, dd, extra_ridge) {
    Ll = lapply(1:ncol(Ba), function(x) Ba[, x])
    Pl = lapply(Ll, function(x) {
        as.spam(t(D) %*% diag(x) %*% D + extra_ridge * diag(1/ncol(Ba), ncol(D)))
    })
    return(list(Pl = Pl, Ll = Ll))
}


#' IWLS step
#' 
#' Update splines coefs via IWLS 
#' @param coefs: current spline coefficients
#' @param Fx: observed residual magnetization (response variable)
#' @param B: B-spline matrix
#' @param Ci: integration matrix
#' @param Pl: adaptive penalty components
#' @param la: current smoothing parameters coefficients
iwls_step = function(coefs, Fx, B, Ci, Pl, la) {
    eta = c(coefs)
    phi = exp(B %*% eta)
    h = phi
    mu = (Ci %*% (h))
    res = c(Fx - mu)

    # Gradient
    H = diag.spam(c(h))

    # set LHS
    Xb = (Ci) %*% H %*% B
    tXb = t(Xb)
    tXbXb = tXb %*% Xb

    # Update pars
    Pprod = lapply(1:length(Pl), function(i) Pl[[i]] * (la[i]))
    P = Reduce("+", Pprod)
    Gi = tXbXb + P
    parmNew = solve(Gi, tXb %*% (res + Xb %*% eta))

    # Ouput
    return(list(par = parmNew, tXbXb = tXbXb, res = res, h = h))
}

#' Update variance components
#' 
#' Update components of the adapitve penalty
#' @param tXbXb: quaddatic form iwls
#' @param eta: current spline coefficients
#' @param res: current working resiquals
#' @param Pl: adaptive penalty components
#' @param la: current smoothing parameters coefficients
#' @return list of variance components
update_variance_components = function(tXbXb, eta, res, Pl, la) {
    lP = append(list(tXbXb), lapply(1:length(Pl), function(i) Pl[[i]]))
    ADcholC = LMMsolver:::ADchol(lP)
    theta = c(1, la)
    dldet = LMMsolver:::dlogdet(ADcholC, theta)
    gk = dldet[-1] * la
    ed = (theta * dldet)[1]

    lP = append(list(diag.spam(ncol(tXbXb))), lapply(1:length(Pl), function(i) Pl[[i]]))
    ADcholC = LMMsolver:::ADchol(lP)
    theta = c(0, la)
    fk = LMMsolver:::dlogdet(ADcholC, theta)[-1] * la

    psi = sapply(Pl, function(x) c(t(eta) %*% x %*% eta))/(fk - gk)
    vphi = sum((res)^2)/(n - ed)
    lan = pmax(1e-08, vphi/(1e-10 + psi))

    return(list(lan = lan, ed = ed, edk = (fk - gk)))
}

#' Adaptive penalty 
#'
#' Create adaptive penalty as sum of its components
#' @param la: vector of smoothing patametes
#' @param Pl: list of penalty components
#' @return adaptive penalty matrix
adapt_pen = function(la, Pl) {
    Pprod = lapply(1:length(Pl), function(i) Pl[[i]] * (la[i]))
    P = Reduce("+", Pprod)
    return(P)
}