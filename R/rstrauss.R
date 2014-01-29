#' Simulate Strauss process in 2- or 3-dimensional box
#' 
#' Simulate the spatial point Strauss process in 2- or 3-dimensional box 
#' (rectangular cuboid). Three algorithms are available:
#' \describe{
#' \item{BD}{Birth-and-death simulation with variable number of points.}
#' \item{dCFTP}{Dominated Coupling-From-The-Past with variable number of points.}
#' \item{MH}{Fixed number of points Metropolis-Hastings.}
#' }
#' 
#' @param beta The beta parameter, controls (but does not equal) intensity
#' @param gamma The repulsion parameter, from 0 to 1 (incl.)
#' @param range Range of pairwise interaction
#' @param n The number of points to distribute, overrides \code{beta}
#' @param bbox A 2 x d  matrix with 1st row giving the lower and 2nd row the higher coordinate ranges for the simulation box.
#' @param iter Number of iterations (in MH & BD) and max. number of steps for dCFTP.
#' @param toroidal Whether to use toroidal correction.
#' @param verb Control verbosity
#' @param perfect Use dCFTP simulation? Otherwise BD, unless n given which use MH.
#' 
#' @details
#' The density of a realisation x of Strauss(beta, gamma, r), where r is the range, is
#' \deqn{f(x)= alfa beta^n(x) gamma^s(x;r)}
#' with scaling constant \code{alfa}, number of points \code{n(x)}, and the number of r-close pairs \code{s(x;r)}.
#' 
#' Under the condition \code{n(x)=n}, the \code{beta==1}.
#' @return
#' A list with
#' \item{x}{2- or 3-dimensional table with the outcome point coordinates}
#' \item{bbox}{the bounding box used for simulation.}
#' 
#' @examples
#' bbox2d <- cbind(c(0,5), c(0, 2.5))
#' bbox3d <- cbind(bbox2d, c(0,2.5))
#' x2 <- rstrauss(beta=10, gamma=0.01, range=0.1, iter=1e5, bbox=bbox2d)
#' x2p <- rstrauss(beta=10, gamma=0.01, range=0.1, perfect=TRUE, bbox=bbox2d, verb=T)
#' x3 <- rstrauss(beta=100, gamma=0.2, range=0.2, perfect=TRUE, bbox=bbox3d, iter=5e4)
#' @import Rcpp
#' @export
#' @useDynLib rstrauss

rstrauss <- function(beta=100,
                     gamma=1, 
                     range=0.1,
                     n,
                     bbox=cbind(c(0,1), c(0,1), c(0,1)), 
                     iter = 1e4,
                     toroidal=FALSE,
                     verb=FALSE,
                     perfect=FALSE) {
  d <- ncol(bbox)
  win <- unlist(bbox)
  # choose algorithm
  # conditional MH simulation
  if(!missing(n)) 
    xyz <- rstrauss_MH(n, gamma, range, win, toroidal, iter, verb)
  # else we have a non-fixed number of points
  else if(perfect) 
    xyz <- rstrauss_DCFTP(beta, gamma, range, win, toroidal, T0=2, dbg=verb, maxtry=iter)
  else 
    xyz <- rstrauss_BD(beta, gamma, range, win, toroidal, iter, verb)  
  # 
  
  xyz <- do.call(cbind, xyz)
  list(x=xyz, bbox=bbox)
}
