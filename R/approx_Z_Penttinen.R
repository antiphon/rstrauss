#' approximate Strauss constant using Penttinen84 approx.
#' 
#' @export 

approximate_strauss_constant_penttinen <- function(beta, gamma, range, bbox, Nmax=10000, n, ...){
  N <- 0:Nmax
  V <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  # for a fixed n:
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  if(missing(n)){
    # Expectation over n:  
    pn <- N*log(beta) - (0.5*N*(N-1)*a/V)
    e <- sum(   exp( -V + N*log(V) +  pn - lfactorial(N)    )  )
    if(!is.finite(e)){
      warning(paste("non-finite value, reducing Nmax to ", Nmax <- round(Nmax * 0.75)) )
      e <- approximate_strauss_constant_penttinen(beta, gamma, range, bbox, Nmax)
    }
    v<- log(e)
  }
  else v <- n*log(V)- (0.5*n*(n-1)*a/V)
  v
}



