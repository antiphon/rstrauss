#' Approximate Strauss constant (master function)
#' 
#' Approximate the (log of) normalizing constant in the Strauss density. Log-scale.
#' 
#' @param beta beta
#' @param gamma gamma
#' @param range range
#' @param bbox bounding box
#' @param method OT, PS, P, R or R3. Default: R3
#' @param ... further details for methods
#' 
#' For all methods but "PS" conditioning on number of points is by setting the named parameter n.
#' 
#' Methods: 
#' 
#' OT : Direct approx., Ogata&Tanemura 1981. \code{\link{approximate_strauss_constant_OT}}
#' 
#' P: Direct approx., Penttinen 1998. \code{\link{approximate_strauss_constant_penttinen}}
#' 
#' PS: MC unbiased estimate, Path sampling Berthelsen&Moller 2003. \code{\link{approximate_strauss_constant_PS}}
#' 
#' R: My own version of Penttinen for unconditional case. Based on a normal approximation (not published). Works when beta*Area >> 0
#' 
#' R3: Third order approx. Same as R but with extra term. (See Ripley 88 p.60-62).\code{\link{approximate_strauss_constant_R3}}
#' 
#' The returned value is in log-scale.
#' 
#' @examples
#' beta <- 100
#' gamma <- 0.1
#' R <- 0.05
#' bbox <- cbind(0:1, 0:1)
#' for(m in c("P", "PS", "R3")) print(approximate_strauss_constant(beta, gamma, R, bbox, m))
#' 
#' approximate_strauss_constant(beta, gamma, R, bbox, n=100)
#' @export

approximate_strauss_constant <- function(beta, gamma, range, bbox, method="R", ...) {
  if(missing(beta)|missing(gamma)|missing(range))stop("Provide beta, gamma and range.")
  if(missing(bbox)) stop("Provide bounding box, column-wise ranges matrix.")
  methods <- c("OT", "P", "PS", "R", "R3")
  fs <- sapply(paste0("approximate_strauss_constant_", c("OT", "penttinen", "PS", "R", "R3") ), get)
  i <- which(match.arg(method, methods)==methods)
  fs[[i]](beta, gamma, range, bbox, ...)
}




