#' approximate Strauss constant using Penttinen84 approx.
#' 
#' @export 

approximate_strauss_constant_penttinen <- function(beta, gamma, range, bbox, Nmax=10000){
  N <- 0:Nmax
  V <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  # for a fixed n:
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  pn <- N*log(beta) - (0.5*N*(N-1)*a/V)
  
  e <- sum(   exp( -V + N*log(V) +  pn - lfactorial(N)    )  )
  log(e)
  #(beta-1)*V  - (beta^2*V*a/2)
}



