#' Approximate Strauss constant (master function)
#' 
#' Approximate the normalizing constant (Z) in the Strauss density.
#' 
#' @param beta beta
#' @param gamma gamma
#' @param range range
#' @param bbox bounding box
#' @param method OT, PS or P
#' @param ... further details for methods
#' 
#' Methods: 
#' OT : Direct approx., Ogata&Tanemura 1981. \code{\link{approximate_strauss_constant_OT}}
#' 
#' P: Direct approx., Penttinen 1998. \code{\link{approximate_strauss_constant_penttinen}}
#' 
#' PS: MC unbiased estimate, Path sampling Berthelsen&Moller 2003. \code{\link{approximate_strauss_constant_PS}}
#' 
#' @export

approximate_strauss_constant <- function(beta, gamma, range, bbox, method="P", ...) {
  if(missing(beta)|missing(gamma)|missing(range))stop("Provide beta, gamma and range.")
  if(missing(bbox)) stop("Provide bounding box, column-wise ranges matrix.")
  methods <- c("OT", "P", "PS")
  fs <- sapply(paste0("approximate_strauss_constant_", c("OT", "penttinen", "PS") ), get)
  i <- which(match.arg(method, methods)==methods)
  fs[[i]](beta,gamma, range, bbox, ...)
}