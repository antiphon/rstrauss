#' approximate unconditional Strauss constant using Gaussian approx. of the expectation
#' 
#' Same Poisson approximation for Z_n as in Penttinen, 
#' but the unconditional is via Gaussian expectation
#' instead of numerical integration.
#' 
#' @export 

approximate_strauss_constant_R <- function(beta, gamma, range, bbox, n){
  V <- prod(apply(bbox, 2, diff))
  dim <- ncol(bbox)
  a <- pi * (1-gamma) * (if(dim==2) (range^2) else ((4/3)*range^3))
  # for a free n:
  if(missing(n)){
    a <- -0.5*a/V
    b <- V*beta
    cee <- 0.5/b-a
    d <- 0.5*(a-1)/cee
    e <- 0.5/cee
    v <- -V + b + 0.5*log(e/b) - 0.5*b + d^2*cee
  } # for a fixed n
  else v <- n*log(V)- (0.5*n*(n-1)*a/V)
  v
}



