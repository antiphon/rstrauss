#' Theoretical maximum for the range parameter for hard core
#' 
#' Given the number of points (intensity) and the gamma=0, what is the 
#' theoretical upper limit for range?
#' 
#' @export

strauss_theory <- function(intensity=100, dim=3) {
  if(dim==3) {
    maxr <- 2*(3/(2*pi^2*intensity))^(1/3)
    c(max_range=maxr)
  }
  else{
    maxr <-2*(2/(pi^2*intensity))^(1/2)
    c(max_range=maxr)
  }
}