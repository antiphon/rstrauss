#' approximate Strauss constant using Ogata&Tanemura approx.
#' 
#' @export 

approximate_strauss_constant_OT <- function(beta, gamma, range, bbox, Nmax=10000){
  N <- 0:Nmax
  winVol <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  pn <- N*log(beta*winVol) + (N*(N-1)/2)*log(1-a/winVol)  - lfactorial(N)
  log(sum(exp(pn)))
}



