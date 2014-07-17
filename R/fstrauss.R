#' Fit stationary and isotropic Strauss process in 2D and 3D
#' 
#' NOTE: This is an early version.
#' 
#' @param x Point pattern, see details for format.
#' @param R The fixed interaction radius.
#' @param approx Which approximation method to use. Default: "LL". See details.
#' @param ... Passed on to the method.
#' 
#' Data can be given in several formats:
#' \describe{
#' \item{matrix}{Columnwise coordinates. Bounding box computed.}
#' \item{ppp}{2D or 3D pattern from spatstat.}
#' \item{list}{Should have elements $x columnwise matrix of coordinates and $bbox bounding box (col-wise ranges).}
#' } 
#' 
#' Three approximation methods are available:
#' \describe{
#' \item{LL}{Logistic likelihood approximation (Baddeley et al 2014). See \link{fstrauss.logistic}}
#' \item{D}{Direct optimization using approximation of the normalizing constant. See \link{fstrauss.direct}}
#' \item{B}{Bayesian version of D. \link{fstrauss.bayes}}
#' }
#' 
#' @export

fstrauss <- function(x, R, approx="D", ...) {
  m <- pmatch(approx, avail<-c("LL","D","B", "VB"))
  if(is.na(m))stop(paste("Method should be one of:", paste(avail, collapse=", ")))
  X <- convert_to_pp(x)
  f <- get(paste0("fstrauss.", c("logistic", "direct", "bayes", "vb"))[m])
  result <- f(x=X, R=R, ...)
  
  #' class
  result$method <- avail[m]
  class(result) <- c("straussfit", class(result))
  result
}

