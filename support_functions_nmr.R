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
#' @param extra_ridge: float for identifiability constraint 
#' @return list of adaptive penalty components
define_adaptive_penalty = function(Ba, D, dd, extra_ridge)
{
    Ll = lapply(1:ncol(Ba), function(x) Ba[, x])
    Pl = lapply(Ll, function(x) {
    as.spam(t(D) %*% diag(x) %*% D + extra_ridge * diag(1 / ncol(Ba), nb))
    })
    return(list(Pl = Pl, Ll = Ll))
}


#' Fit function
#' 
#' Function to fit smooth deconvolution model
#' @param Fx: observed residual magnetization
#' @param Ci: deconvolution matrix for integration
#' @param Pl: adaptive penalty components
#' @param Ll: list with a columns of the adaptive penalty bases
#' @param maxkit: max number of updates for smoothing parameter estimation
#' @param maxit: maximum number of internal IWLS steps (coefficient estimation)
#' @return list of results #TODO
fit_smooth_deconvolution = function(Fx, Ci, Pl, Ll, B, maxkit = 200, maxit = 100){    
    # Initialize
    nci = nrow(B)
    nb  = ncol(B)
    la  = 10^rep(2, ncol(Ba))
    eta = solve(t(Ci %*% B) %*% (Ci %*% B) + 1e2 * t(D) %*% D, t(Ci %*% B) %*% Fx, tol = 1e-50
)
    lz  = 0.75
    tol = 1e-5
    param = eta

    for(kit in 1:maxkit)
    {
        # Estimate eta/lambdas
        for(it in 1:maxit)  
        {
            # Set RHS
            eta = c(param)
            phi = exp(B %*% eta) 
            h   = phi 
            mu  = (Ci %*% (h)) 
            res = c(Fx - mu)
        
            # Gradient
            H    = diag.spam(c(h))

            # set LHS
            Xb   = (Ci) %*% H %*% B
            tXb  = t(Xb) 
            tXbXb = tXb %*% Xb
        
            # Update pars
            Pprod = lapply(1:length(Ll), function(i) Pl[[i]]* (la[i]) )
            P  = Reduce('+', Pprod)
            Gi = tXbXb + P 
            parmNew = solve(Gi, tXb %*% (res + Xb %*% eta))
        
            # Check conv
            de    = sum((parmNew - param)**2) / (.1 + sum(param**2))
            print(de)
            crit  = de < tol
            parmNew = lz * parmNew + (1-lz) * param
            if(crit) break
            param   = parmNew
            it      = it + 1
        }    

        # Update variance components
        lP = append(list(tXbXb), lapply(1:length(Ll), function(i) Pl[[i]])) 
        ADcholC = LMMsolver:::ADchol(lP)
        theta = c(1.0, la)
        dldet = LMMsolver:::dlogdet(ADcholC, theta) 
        gk = dldet[-1] * la
        ed = (theta * dldet)[1]

        lP = append(list(diag.spam(nb)), lapply(1:length(Ll), function(i) Pl[[i]]))
        ADcholC = LMMsolver:::ADchol(lP)
        theta = c(0, la)
        fk = LMMsolver:::dlogdet(ADcholC, theta)[-1] * la

        psi = sapply(Pl, function(x) c(t(eta) %*% x %*% eta)) / (fk - gk)
        vphi= sum((res)**2) / (n - ed) 
        lan = pmax(1e-8, vphi/(1e-10 + psi))

        # Check convergence
        crit_kp = sum((la-lan)**2) / (.1 + sum(la**2))
        la      = lan 
        cat("iter Outer:", kit, "crit:", crit_kp, "iter Inner:", it, "\n")

        if(crit_kp < tol) break 
    }
    
    
    out = list(eta = eta, h = h, res = res, la = la, ed = ed, edk = (fk - gk), kit = kit, conv = kit < maxkit)
    return(out)
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
