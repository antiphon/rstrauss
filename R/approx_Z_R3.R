#' approximate Strauss constant using Gaussian approx. of the expectation
#' 
#' Conditional and unconditional approximation of normalizing constant
#' of the stationary homogeneous Strauss model.
#' Third order approximation, as Penttinen is second order.
#' 
#' @param beta Beta
#' @param gamma Gamma
#' @param range Interaction range
#' @param bbox Bounding box (col. wise ranges), used for dimension and volume.
#' @param n if given, compute the conditional Z.
#' @param ... ignored.
#' 
#' Expectation to get the unconditional is using Gaussian approximation.
#' 
#' This is the best in the package.
#' 
#' @export 

approximate_strauss_constant_R3 <- function(beta, gamma, range, bbox, n, ...){
  V <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  #' dimension dependent
  Cr <- pi * ifelse(dim==2, range^2, 4*range^3/3)
  c0 <- ifelse(dim==2, pi * (4*pi-3*sqrt(3))/4, 5*pi^3/6) 
  #' constant
  a1 <- log(V)
  a2 <- (gamma-1) * Cr/(2*V)
  a3 <- c0 * (gamma-1)^3 * range^(dim*2)/(6*V^2)
  #'
  #' for a free n:
  if(missing(n)){
    k <- V*beta
    kk <- 1/(1/k - 2*a2)
    mu <- (1+a1)*kk
    v <- V*(beta-1) + mu^2/(2*kk) - k/2 + 0.5*log(kk/k) + a3*(mu^3 + 3*mu*kk )
  } # for a fixed n
  else{
    v <- a1 * n + a2 * n^2 + a3 * n^3
  } 
  v
}

