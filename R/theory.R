#' Theoretical maximum for the range parameter for hard core
#' 
#' Given the number of points (intensity) and the gamma=0, what is the 
#' theoretical upper limit for range?
#' 
#' @export

strauss_theory <- function(intensity, dim=3) {
  if(dim==3) {
    maxr <- 2*(3/(2*pi^2*intensity))^(1/3)
    c(max_range=maxr)
  }
  else{
    maxr <-2*(2/(pi^2*intensity))^(1/2)
    c(max_range=maxr)
  }
}

#' Approximate intensity for Strauss based on beta,gamma,range
#' 

#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param R range
#' @param dim dimension
#' @param stucki see details
#'
#' @details 
#' if stucki=FALSE, use Poisson mgf approximation, 
#' else return the limits given in Stucki 2012 paper, which are 
#' accurate but bounds. The average is quite good.
#' 
#' @export

strauss_intensity <- function(beta, gamma, range, dim=3, stucki=FALSE){
  G <- (1-gamma) * range^dim * ( if(dim==3) 3*pi/4 else pi )
  if(!stucki){
    LambertW(beta*G)/G
  }else{
    c(lo=beta/(1+beta*G), hi=beta/(2-exp(-beta*G)))
  }
}

#' Solve beta for Strauss given on n,gamma,range
#' 
#' @param n number of points 
#' @param gamma gamma parameter
#' @param R range
#' @param dim dimension
#' @param stucki see details
#' @param take see details
#' 
#' @details
#' If stucki=FALSE, use Poisson mgf approximation, 
#' else use the function 'take' (default mean) on the bounds given by Stucki 2012.
#' 'take' could be mean, max or min, for example.
#' see \link{strauss_intensity}
#' 
#' @export

strauss_beta <- function(n, gamma, R, dim=3, stucki=FALSE, take=mean){
  f <- if(stucki) function(x) (strauss_intensity(x, gamma, R, dim) - n)^2
        else function(x) (take(strauss_intensity(x, gamma, R, dim, TRUE)) - n)^2
  nlm(f, n)$estimate
}
