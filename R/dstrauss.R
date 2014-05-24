#' density of a strauss process
#' 
#' Now only for cuboidal window
#' 
#' @param x list with x: coordinate matrix (col-wise), bbox: bounding box (col-wise ranges)
#' @param beta beta parameter
#' @param gamma gamma parameter
#' @param range interaction range
#' @param method how to approximate the constant.
#' @param ... passed to constant approximator.
#' 
#' See \code{\link{approx_strauss_constant}} for the methods
#' 
#' @import cutgeom
#' @export

dstrauss <- function(x, beta, gamma, range, method="OT", ...) {
  pairs <- sum( sapply(geom(x$x, r=range), length)  )/2
  winV <- prod(apply(x$bbox, 2, diff))
  dim <- ncol(x$bbox)
  #' the constant:
  Z <- approximate_strauss_constant_OT(beta, gamma, range, bbox=x$bbox, method=method, ...)
  #'
  beta^nrow(x$x) * gamma^pairs / Z
}

