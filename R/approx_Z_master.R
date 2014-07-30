#' Approximate Strauss constant (master function)
#' 
#' Approximate the normalizing constant (Z) in the Strauss density.
#' 
#' @param beta beta
#' @param gamma gamma
#' @param range range
#' @param bbox bounding box
#' @param method OT, PS, P or R. Default: R
#' @param ... further details for methods
#' 
#' Conditioning on number of points can be done with OT and P by setting the named parameter n.
#' 
#' Methods: 
#' 
#' OT : Direct approx., Ogata&Tanemura 1981. \code{\link{approximate_strauss_constant_OT}}
#' 
#' P: Direct approx., Penttinen 1998. \code{\link{approximate_strauss_constant_penttinen}}
#' 
#' PS: MC unbiased estimate, Path sampling Berthelsen&Moller 2003. \code{\link{approximate_strauss_constant_PS}}
#' 
#' R: My own version for unconditional Z. Based on a normal approximation (not published). Works when beta*Area >> 0
#' 
#' @export

approximate_strauss_constant <- function(beta, gamma, range, bbox, method="R", ...) {
  if(missing(beta)|missing(gamma)|missing(range))stop("Provide beta, gamma and range.")
  if(missing(bbox)) stop("Provide bounding box, column-wise ranges matrix.")
  methods <- c("OT", "P", "PS", "R")
  fs <- sapply(paste0("approximate_strauss_constant_", c("OT", "penttinen", "PS", "R") ), get)
  i <- which(match.arg(method, methods)==methods)
  fs[[i]](beta,gamma, range, bbox, ...)
}