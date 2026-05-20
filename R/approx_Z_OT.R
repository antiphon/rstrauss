#' approximate Strauss constant using Ogata&Tanemura approx.
#' 
#' @export 

approximate_strauss_constant_OT <- function(beta, gamma, range, bbox, Nmax=10000, n, deg=2, ...){
  N <- 0:Nmax
  V <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  if(deg==3){
    
  }
  if(missing(n)){
    pn <- N*(log(beta)+log(V)) + (N*(N-1)/2)*log(1-a/V)  - lfactorial(N)
    v <- log(sum(exp(pn))) - V
    if(!is.finite(v)){
      warning(paste("non-finite value, reducing Nmax to ", Nmax <- round(Nmax * 0.75)) )
      v <- approximate_strauss_constant_OT(beta, gamma, range, bbox, Nmax, n, deg)
    }
  }
  else v <- n*log(V)+(n*(n-1)/2)*log(1-a/V)
  v
}



