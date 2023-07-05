# Install/load libraries
libs = c('colorout', 'parallel', 'R.matlab', 'statmod', 'sfsmisc', 'fields', 'JOPS', 'matrixStats', 'scales', 'ggplot2') 
install.packages(setdiff(libs, rownames(installed.packages())))  
pkgLoad = lapply(libs, function( x ) require(x, character.only = T))    

#' Vector of T2 times
#' 
#' Creates a fector of log10 T2 values
#' @param m: integer length of the T2 times
#' @param min_t: min T2 on log10 scale
#' @param max_t: max T2 on log10 scale
#' @return Vector of T2 values
define_t2_vector = function(m, min_t, max_t)
{
    tt = seq(min_t, max_t, len = m)
    return(tt)
}

#' Model matrix
#' 
#' Computes the model matrix for integral approx.
#' @param x: vector of boserved times (dim n)
#' @param tt: vector of m T2 times (log10 scale)
#' @return An integration matrix (n x m)
define_model_matrix = function(x, tt)
{
    Ci = exp(outer(-x, 1 / (10^tt)))
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
define_adaptive_penalty = function(Ba, D, dd, extra_ridge)
{
    Ll = lapply(1:ncol(Ba), function(x) Ba[, x])
    Pl = lapply(Ll, function(x) {
    as.spam(t(D) %*% diag(x) %*% D + extra_ridge * diag(1 / ncol(Ba), ncol(D)))
    })
    return(list(Pl = Pl, Ll = Ll))
}

#' Fit initialization
#' 
#' Function to initialize fit
#' @param Fx: observed residual magnetization
#' @param nla: number of bases for the adaptive penalty
#' @param Ci: deconvolution matrix for integration
#' @param B: B-spline matrix
#' @param D: Difference penalty matrix
#' @return initial smoothing and spline coefficients
initialize_smooth_deconvolution = function(Fx, nla, Ci, B, D) 
{
    la  = 10^rep(2, nla)
    eta = solve(t(Ci %*% B) %*% (Ci %*% B) + 1e2 * t(D) %*% D, t(Ci %*% B) %*% Fx, tol = 1e-50)
    return(list(la = la, eta = eta))
}

#' Fit function
#' 
#' Function to fit smooth deconvolution model
#' @param Fx: observed residual magnetization
#' @param Ci: deconvolution matrix for integration
#' @param Pl: adaptive penalty components
#' @param B: B-spline matrix
#' @param initials: list of initials spline and smoothing coefficients
#' @param maxkit: max number of updates for smoothing parameter estimation
#' @param maxit: maximum number of internal IWLS steps (coefficient estimation)
#' @param lz: spline coefficients update dumping parameter to improve convergence
#' @param tol: tolerance parameter
#' @return list of results.
fit_smooth_deconvolution = function(Fx, Ci, Pl, B, initials, maxkit = 200, maxit = 100, lz  = 0.75, tol = 1e-5){    
    # Initialize
    nci = nrow(B)
    nb  = ncol(B)
    param = initials$eta
    la = initials$la
    for(kit in 1:maxkit)
    {
        # Estimate eta/lambdas
        for(it in 1:maxit)  
        {
            # Set RHS
            iwls_est = iwls_step(param, Fx, B, Ci, Pl, la)
            parmNew = iwls_est$par
            tXbXb = iwls_est$tXbXb
            res = iwls_est$res

            # Check conv
            de  = sum((parmNew - param)**2) / (.1 + sum(param**2))
            print(de)
            crit  = de < tol
            parmNew = lz * parmNew + (1-lz) * param
            if(crit) break
            param   = parmNew
            it      = it + 1
        }    
        # Update smoothing parameters
        varicance_updates = update_variance_components(tXbXb, param, res, Pl, la)
        lan = varicance_updates$lan

        # Check convergence
        crit_kp = sum((la-lan)**2) / (.1 + sum(la**2))
        la = lan 
        cat("iter Outer:", kit, "crit:", crit_kp, "iter Inner:", it, "\n")

        if(crit_kp < tol) break 
    }
    
    # Output
    h = iwls_est$h
    ed = varicance_updates$ed
    edk = varicance_updates$edk
    out = list(eta = param, h = h, res = res, la = la, ed = ed, edk = edk, kit = kit, conv = kit < maxkit)
    return(out)
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
iwls_step = function(coefs, Fx, B, Ci, Pl, la)
{
    eta = c(coefs)
    phi = exp(B %*% eta) 
    h   = phi 
    mu  = (Ci %*% (h)) 
    res = c(Fx - mu)
        
    # Gradient
    H = diag.spam(c(h))

    # set LHS
    Xb = (Ci) %*% H %*% B
    tXb = t(Xb) 
    tXbXb = tXb %*% Xb
        
    # Update pars
    Pprod = lapply(1:length(Pl), function(i) Pl[[i]]* (la[i]))
    P  = Reduce('+', Pprod)
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
update_variance_components = function(tXbXb, eta, res, Pl, la)
{
    lP = append(list(tXbXb), lapply(1:length(Pl), function(i) Pl[[i]])) 
    ADcholC = LMMsolver:::ADchol(lP)
    theta = c(1.0, la)
    dldet = LMMsolver:::dlogdet(ADcholC, theta) 
    gk = dldet[-1] * la
    ed = (theta * dldet)[1]

    lP = append(list(diag.spam(ncol(tXbXb))), lapply(1:length(Pl), function(i) Pl[[i]]))
    ADcholC = LMMsolver:::ADchol(lP)
    theta = c(0, la)
    fk = LMMsolver:::dlogdet(ADcholC, theta)[-1] * la

    psi = sapply(Pl, function(x) c(t(eta) %*% x %*% eta)) / (fk - gk)
    vphi= sum((res)**2) / (n - ed) 
    lan = pmax(1e-8, vphi/(1e-10 + psi))

    return(list(lan = lan, ed = ed, edk = (fk - gk)))
}

#' Adaptive penalty 
#'
#' Create adaptive penalty as sum of its components
#' @param la: vector of smoothing patametes
#' @param Pl: list of penalty components
#' @return adaptive penalty matrix
adapt_pen = function(la, Pl)
{
    Pprod = lapply(1:length(Pl), function(i) Pl[[i]] * (la[i]) )
    P = Reduce('+', Pprod)
    return(P)
}

#' Smooth deconvolution 
#' 
#' Fit a smooth deconvolution model to observed relaxation
#' @param Fx: observed residual magnetization (response variable)
#' @param bdeg: B-spline degree
#' @param ndx: number of internal B-spline knots
#' @param dd: order of the difference penalty
#' @param min_t: minimum u-value (log10 scale) for integration grid
#' @param max_t: maximum u-value (log10 scale) for integration grid
#' @param m: number of intrgration point
#' @param extra_ridge: extra penalty to condition P
#' @return list of results.
smooth_deconvolution = function(Fx, bdeg = 3, ndx = 35, dd = 2, min_t = -3.5, max_t = 3.5, m = 200, extra_ridge=1e-12)
{
    # Define model matrix (with integration rule)
    tt = define_t2_vector(m, min_t = min_t, max_t = max_t)
    Ci = define_model_matrix(x, tt)

    # Define bases and penalty
    B = as.spam(bbase(tt, bdeg = bdeg, nseg = ndx))
    nb = ncol(B)
    D = diff(diag(nb), diff = dd)

    # Bases for adaptive penalty
    xa = 1:(nb - dd) / (nb - dd)
    Ba = bbase(xa, nseg = 5, bdeg = 3)
    pen_comp = define_adaptive_penalty(Ba, D, dd, extra_ridge=extra_ridge)
    Pl = pen_comp$Pl

    # Initialization
    scl = max(Fx)
    initials = initialize_smooth_deconvolution(Fx/scl, ncol(Ba), Ci, B, D)

    # Estimate eta and lambda
    ests = fit_smooth_deconvolution(Fx = Fx/scl, Ci = Ci, Pl = Pl, B = B, initials = initials)
    h = ests$h * scl
    la = ests$la
    mu = Ci %*% h 
    res = ests$res * scl

    return(list(h = h, la=la, lla = log10(Ba %*% la), res = res, mu = mu, tt = tt, xa = xa))
}
