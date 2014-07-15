#' approximate Strauss constant using Path sampling
#' 
#' @param nsim number of sims for MC expectation
#' @param steps grid reso for gamma
#' @param perfMaxIter max iter for the first perfect pattern
#' @param MHiter iterations per MH simulations starting from perfect pattern
#' @param ... passed on to rstrauss simulator
#' @return
#' Value is the natural log of Z(theta)
#' @import multicore cutgeom
#' @export 

approximate_strauss_constant_PS <- function(beta, gamma, range, bbox, nsim=10, steps=10, perfMaxIter=1e5, MHiter=1000, ...){
  if(steps<3) stop("Need gamma grid steps > 2")
  #'
  #' Simulate in a dilated bbox to encounter edge effects
  bbox_sim <- bbox_dilate(bbox, range)
  #'
  #'helpers: Simulate a pattern given gamma
  rst0 <- function(gam) rstrauss(beta=beta, gamma=gam, range=range, bbox=bbox_sim, perfect=TRUE, iter=perfMaxIter, ...)
  rst <- function(x0, gam) rstrauss(n=nrow(x0$x), gamma=gam, range=range, bbox=x0$bbox, start=x0$x, iter=MHiter, ...)
  #'
  #' TODO: fix the border correction
  #' Calcuate one Sr, use reduced window border correction
  Sr <- function(x) {
    inside <- which( bbox_distance(x$x, x$bbox) > range )
    outside <- setdiff(1:nrow(x$x), inside)
    G <- geom(x$x, from=inside, to=inside, r=range)
    G2 <- geom(x$x, from=inside, to=outside, r=range)
    #'
    
    sum(sapply(G, length))/2 + sum(sapply(G2, length)) # number of pairs
  }
  #'
  #' path integral grid
  ggrid <- seq(gamma, 1, length=steps)
  dg <- diff(ggrid[1:2])
  #' The expectations
  E <- NULL
  dummy <- as.list(2:nsim)
  #' loop for estimating the expectations
  for(g in ggrid){
    sim0 <- rst0(g)
    sf <- function(v) rst(sim0, g)
    sims <- mclapply(dummy, sf)
    Srs <- c(Sr(sim0), sapply(sims, Sr))
    Ex <- sum(Srs)/nsim
    E <- c(E, Ex)
  }
  #' then approximate the integral
  bw <- 2:(steps-1) # intermediate values on the path
  A <- dg * ( E[1]/(2*gamma) + E[steps]/2 + sum(E[bw]/ggrid[bw]) )
  #' then back to Z, log scale
  V <- prod(apply(bbox, 2, diff))
  (beta-1)*V - A
}

