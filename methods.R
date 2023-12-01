source("utils.R")

#' Initialization
#' 
#' Initialize class smooth_deconv
setup_smooth_deconv = function(Fx, x, m, min_t, max_t, dd, ndx, bdeg, extra_ridge) {
    obj = list(Fx = Fx, scl = max(Fx), x = x, m = m, min_t = min_t, max_t = max_t,
        dd = dd, ndx = ndx, bdeg = bdeg, extra_ridge = extra_ridge)
    class(obj) = "smooth_deconv"
    return(obj)
}


#' Fit initialization
#' 
#' Function to initialize fit
#' @param obj: object of class smooth_deconv'
#' @return initial smoothing and spline coefficients
initialize.smooth_deconv = function(obj) {
    m = obj$m
    min_t = obj$min_t
    max_t = obj$max_t
    bdeg = obj$bdeg
    ndx = obj$ndx
    dd = obj$dd
    extra_ridge = obj$extra_ridge

    # Define model matrix (with integration rule)
    tt = define_t2_vector(m, min_t, max_t)
    Ci = define_model_matrix(x, tt)

    # Define bases and penalty
    B = spam::as.spam(JOPS::bbase(tt, bdeg = bdeg, nseg = ndx))
    nb = ncol(B)
    D = diff(diag(nb), diff = dd)

    # Bases for adaptive penalty
    xa = 1:(nb - dd)/(nb - dd)
    Ba = JOPS::bbase(xa, nseg = 5, bdeg = 3)
    pen_comp = define_adaptive_penalty(Ba, D, dd, extra_ridge = extra_ridge)
    Pl = pen_comp$Pl

    # Coefficients
    nla = ncol(Ba)
    la = 10^rep(2, nla)
    scl = obj$scl
    eta = solve(t(Ci %*% B) %*% (Ci %*% B) + 100 * t(D) %*% D, t(Ci %*% B) %*% (Fx/scl),
        tol = 1e-50)

    # Output
    obj$la = la
    obj$eta = eta
    obj$Ci = Ci
    obj$Pl = Pl
    obj$B = B
    obj$Ba = Ba
    obj$tt = tt
    obj$xa = xa

    return(obj)
}

#' Fit function
#' 
#' Function to fit smooth deconvolution model
#' @param obj: an object of class 'smooth_deconv'
#' @param maxkit: max number of updates for smoothing parameter estimation
#' @param maxit: maximum number of internal IWLS steps (coefficient estimation)
#' @param lz: spline coefficients update dumping parameter to improve convergence
#' @param tol: tolerance parameter
#' @return list of results.
fit.smooth_deconv = function(obj, maxkit = 200, maxit = 100, lz = 0.75, tol = 1e-05) {

    # Initialize
    param = obj$eta
    la = obj$la
    Fx = obj$Fx/obj$scl
    B = obj$B
    Ci = obj$Ci
    Pl = obj$Pl
    nci = nrow(B)
    nb = ncol(B)

    # Loop
    for (kit in 1:maxkit) {
        for (it in 1:maxit) {
            iwls_est = iwls_step(param, Fx, B, Ci, Pl, la)
            parmNew = iwls_est$par
            tXbXb = iwls_est$tXbXb
            res = iwls_est$res

            de = sum((parmNew - param)^2)/(0.1 + sum(param^2))
            print(de)
            crit = de < tol
            parmNew = lz * parmNew + (1 - lz) * param
            if (crit)
                break
            param = parmNew
            it = it + 1
        }
        variance_updates = update_variance_components(tXbXb, param, res, Pl, la)
        lan = variance_updates$lan

        crit_kp = sum((la - lan)^2)/(0.1 + sum(la^2))
        la = lan
        cat("iter Outer:", kit, "crit:", crit_kp, "iter Inner:", it, "\n")

        if (crit_kp < tol)
            break
    }

    # Output
    h = iwls_est$h
    ed = variance_updates$ed
    edk = variance_updates$edk
    obj$eta = param
    obj$h = h * obj$scl
    obj$res = res * obj$scl
    obj$la = la
    obj$ed = ed
    obj$edk = edk
    obj$kit = kit
    obj$conv = kit < maxkit
    return(obj)
}

#' Plot smooth deconvolution
#' 
#' Plots the resutls of a smooth deconvoultion fit
#' @param obj: an pbject of class smooth_deconv
#' @return ggplot with 4 panels: data and fit, spectrum,
#' residuals and adaptive penalty
plot.smooth_deconv = function(obj) {
    h = obj$h
    la = obj$la
    resid = obj$resid
    mu = obj$mu
    res = obj$res
    lla = obj$lla
    tt = obj$tt
    xa = obj$xa

    # Plot results
    fitDat = data.frame(x = x, obs = Fx, fit = mu, resid = res)
    estDat = data.frame(time = 10^tt, h = h)
    penDat = data.frame(xa = xa, pen = lla, ta = 10^seq(min(tt), max(tt), len = length(xa)))

    fitPlt = ggplot(fitDat, aes(x = x, y = obs)) + geom_point() + geom_line(aes(x = x,
        y = fit), colour = "orange") + xlab("Time (ms)") + ylab("Residual magnetisation") +
        ggtitle(paste("Sample", ix, "(DM =", paste0(round(data$Y[ix, 1], 5), "%)"))) +
        theme_bw()

    estPlt = ggplot(estDat, aes(x = time, y = h)) + geom_line() + scale_x_continuous(trans = "log10",
        limits = range(data$Time) + c(-0.6, 1300)) + xlab("Time (ms) - log10 scale") +
        ylab("Density") + theme_bw()

    resPlt = ggplot(fitDat, aes(x = x, y = resid)) + geom_point() + geom_hline(yintercept = 0,
        colour = "gray50") + xlab("Time (ms)") + ylab("Residuals") + theme_bw()

    penPlt = ggplot(penDat, aes(x = ta, y = pen)) + geom_point() + geom_line() +
        scale_x_continuous(trans = "log10") + xlab("Time (ms) - log10 scale") + ylab("Adaptive penalty - log10 scale") +
        theme_bw()

    plRes = gridExtra::grid.arrange(fitPlt, estPlt, resPlt, penPlt, nrow = 2, ncol = 2)
    return(plRes)
}

#' Smooth deconvolution 
#' 
#' Fit a smooth deconvolution model to observed relaxation
#' @param Fx: observed residual magnetization (response variable)
#' @param x: observed time vector
#' @param bdeg: B-spline degree
#' @param ndx: number of internal B-spline knots
#' @param dd: order of the difference penalty
#' @param min_t: minimum u-value (log10 scale) for integration grid
#' @param max_t: maximum u-value (log10 scale) for integration grid
#' @param m: number of intrgration point
#' @param extra_ridge: extra penalty to condition P
#' @return list of results.
smooth_deconvolution = function(Fx, x, bdeg = 3, ndx = 35, dd = 2, min_t = -3.5,
    max_t = 3.5, m = 200, extra_ridge = 1e-12) {

    smooth_deconv = setup_smooth_deconv(Fx, x, m, min_t, max_t, dd, ndx, bdeg, extra_ridge)

    # Initialization
    initials = initialize(smooth_deconv)

    # Estimate eta and lambda
    ests = fit(initials)
    ests$mu = ests$Ci %*% ests$h
    ests$lla = log10(ests$Ba %*% ests$la)

    return(ests)
}


