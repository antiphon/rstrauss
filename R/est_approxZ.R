#' Fit Strauss process to data using maximum approx. likelihood
#' 
#' ML estimate of beta and gamma using optimization of log-likelihood,
#' where the normalizing constant is approximated directly.
#' 
#' @param x a list with elements $x:data.frame of locations (ncol 2 or 3), $bbox: matrix with columns giving window bounds
#' @param R the fixed known interaction range. Use profile pl to infer this.
#' @param lower Lower limits for searching (beta,gamma).
#' @param upper Upper limits for searching (beta,gamma).
#' @param init Initial values for optim.
#' @param ... passed on to \code{\link{approximate_strauss_constant}}
#' 
#' Reduced sample border correction with radius R.
#' 
#' See \link{approximate_strauss_constant} for approximation options.
#' 
#' @export

fstrauss.direct <- function(x, R, lower=c(1e-9, 1e-9), upper=c(Inf, 1), init, ...) {
  bbox <- if(is.list(x)) x$bbox else apply(x, 2, range)
  x <- ( if(is.list(x)) x$x else x )
  dim <- ncol(bbox)
  VR <- R^dim * ( if(dim==3) 3*pi/4 else pi )
  #' reduced sample border correction
  bbox_dil <- bbox_erode(bbox, R)
  V <- prod(apply(bbox_dil, 2, diff))
  #' compute statistics, reduced sample border correction
  bbdists <- bbox_distance(x, bbox)
  inside <- which(bbdists > R)
  outside <- which(bbdists < R)
  N <- length(inside)
  nlist_in <- geom(x, from=inside, to=inside, r=R)
  nlist_2out <- geom(x, from=inside, to=outside, r=R)
  degs_in <- sapply(nlist_in, length)
  degs_2out <- sapply(nlist_2out, length)
  tv <- cbind(N, sum(degs_in)/2 + sum(degs_2out))
  #' initial values
  initial_values <- if(missing(init)){
    lambda <-  N/V
    g0 <-  tv[2]/(N*(N-1)*VR/(2*V))
    G <-  (1 - g0) * VR
    c(exp(lambda*G)*lambda, g0)
  } else init
  
  #' optim function
  fun <- function(theta) {
    lz <- approximate_strauss_constant(theta[1], theta[2], R, bbox_dil, ...)
    v<-lz - tv%*%log(theta) # to maximize must be inverted
    if(is.finite(v)) v else 1e9
  }
  #'
  res <- optim(initial_values, fun, lower=lower, upper=upper,  method="L-BFGS-B")
  #'
  coef  <- c(res$par, R)
  names(coef) <- c("beta", "gamma", "r_given")
  #' log-likelihood
  lz <- approximate_strauss_constant(coef[1], coef[2], R, bbox_dil, ... )
  ll <- tv%*%log(coef[1:2])-lz
  list(theta=coef, optim_out=res, logLik=ll)
}


