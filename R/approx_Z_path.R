#' approximate Strauss constant using Path sampling
#' 
#' @param nsim number of sims for MC expectation
#' @param steps grid reso for gamma
#' @param ... passed on to rstrauss simulator
#' @return
#' Value is the natural log of Z(theta)
#' @import multicore cutgeom
#' @export 

approximate_strauss_constant_PS <- function(beta, gamma, range, bbox, nsim=10, steps=10, ...){
  if(steps<3) stop("Need gamma grid steps > 2")
  #'
  #'helpers: Simulate a pattern given gamma
  rst <- function(gam) rstrauss(beta=beta, gamma=gam, range=range, bbox=bbox, ...)
  #'
  #' TODO: fix the border correction
  #' Calcuate one Sr, use reduced window border correction
  Sr <- function(x) {
    G <- geom(x$x, r=range)    
    #ok <- 1:length(G) #bbox_distance(x$x, bbox) > range # exclude too close to border
#     print(sum(ok)/length(ok))
    #sum(sapply(G[ok], length))/2 # number of pairs
    sum(sapply(G, length))/2 # number of pairs
  }
  #' path integral grid
  ggrid <- seq(gamma, 1, length=steps)
  dg <- diff(ggrid[1:2])
  #' The expectations
  E <- NULL
  dummy <- as.list(1:nsim)
  #' loop for estimating the expectations
  for(g in ggrid){
    sf <- function(v) rst(g)
    sims <- mclapply(dummy, sf)
    Srs <- sapply(sims, Sr)
    Ex <- sum(Srs)/nsim
    E <- c(E, Ex)
  }
  #' then approximate the integral
  bw <- 2:(steps-1) # intermediate values on the path
  A <- dg * ( E[1]/(2*gamma) + E[steps]/2 + sum(E[bw]/ggrid[bw]) )
  #' then back to Z, log
  beta - A
}
