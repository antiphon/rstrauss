#' Bayesian inference for stationary Strauss process
#' 
#' Using direct approximation of the normalizing constant.
#' 
#' @param x a list with elements $x:data.frame of locations (ncol 2 or 3), $bbox: matrix with columns giving window bounds
#' @param R the fixed known interaction range. Use profile pl to infer this.
#' @param priors A list of priors. See details.
#' @param lower Lower limits for searching (beta,gamma).
#' @param upper Upper limits for searching (beta,gamma).
#' @param init Initial values for optim.
#' @param ... passed on to \code{\link{approximate_strauss_constant}}
#' 
#' Reduced sample border correction with radius R.
#' 
#' We use Gamma(a,b) prior for beta and Beta(c,d) prior for Gamma (funny...)
#' so priors-parameter should be a list with elements $gamma=c(a,b) and $beta=c(c,d).
#' 
#' At the moment only the MAP is computed and returned in the $par element.
#' 
#' See \link{approximate_strauss_constant} for approximation options.
#' 
#' 
#' @export

fstrauss.bayes <- function(x, R, priors=list(beta=c(1,1e-5), gamma=c(1,1)), 
                           lower=c(1e-9, 1e-9), upper=c(Inf, 1), init, 
                           ...) {
  warning("Bayesian fitter not fully implemented, only MAP estimate available.")
  bbox <- if(is.list(x)) x$bbox else apply(x, 2, range)
  bbox_dil <- bbox_erode(bbox, R)
  V <- prod(apply(bbox_dil, 2, diff))
  x <- ( if(is.list(x)) x$x else x )
  dim <- ncol(bbox)
  VR <- R^dim * ( if(dim==3) 3*pi/4 else pi )
  # compute statistics, reduced sample border correction
  bbdists <- bbox_distance(x, bbox)
  inside <- which(bbdists > R)
  outside <- which(bbdists < R)
  N <- length(inside)
  nlist_in <- geom(x, from=inside, to=inside, r=R)
  nlist_2out <- geom(x, from=inside, to=outside, r=R)
  degs_in <- sapply(nlist_in, length)
  degs_2out <- sapply(nlist_2out, length)
  tv <- cbind(N, sum(degs_in)/2 + sum(degs_2out))
  #
  initial_values <- if(missing(init)){
    lambda <-  N/V
    g0 <-  tv[2]/(N*(N-1)*VR/(2*V))
    G <-  (1 - g0) * VR
    c(exp(lambda*G)*lambda, g0)
  } else init
  # optim function, includes priors
  fun <- function(theta) {
    lz <- approximate_strauss_constant(theta[1], theta[2], R, bbox_dil, ...)
    lt <- log(theta)
    lpri <- lt[1]*(priors$beta[1]-1)-priors$beta[2]*lt[2] + 
            lt[2]*(priors$gamma[1]-1)+priors$gamma[2]*log(1-theta[2])
    
    v<-lz - tv%*%log(theta) - lpri # to maximize must be inverted
    if(is.finite(v)) v else 1e9
  }
  #
  res <- optim(c(N/V, 0.5), fun, lower=lower, upper=upper,  method="L-BFGS-B")
  #
  coef  <- c(res$par, R)
  names(coef) <- c("beta", "gamma", "r_given")
  
  list(theta=coef, optim_out=res, logLik=NA)
}