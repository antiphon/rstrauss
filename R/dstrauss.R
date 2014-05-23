#' density of a strauss process
#' 
#' Now only for cuboidal window
#' 
#' @param x list with x: coordinate matrix (col-wise), bbox: bounding box (col-wise ranges)
#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param range interaction range
#' @param method how to approximate.
#' 
#' Currently only method="OT" for Ogata&Tanemura 1981 is available, using numerical summation. 
#' 
#' @import cutgeom
#' @export

dstrauss <- function(x, beta, gamma, range, method="OT") {
  pairs <- sum( sapply(geom(y$x, r=range), length)  )/2
  #' the constant:
  Z <- approximate_strauss_constant_OT(beta, gamma, current)
  #'
  beta^sum(current$z) * gamma^pairs / Z
}


#' approximate Strauss constant using Ogata&Tanemura approx.
#' 
#' @export 

approximate_strauss_constant_OT <- function(beta, gamma, range, winV, dim, Nmax=10000){
  N <- 1:Nmax
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  sum( winVol^N * (1-a/winVol)^(N*(N-1)/2) )
}



