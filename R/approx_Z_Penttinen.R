#' approximate Strauss constant using Penttinen84 approx.
#' 
#' @export 

approximate_strauss_constant_penttinen <- function(beta, gamma, range, bbox, Nmax=10000){
  N <- 0:Nmax
  winVol <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  pn <- N*log(beta) - (0.5*N*(N-1)*a)  - lfactorial(N)
  log(sum(exp(pn)))
}



