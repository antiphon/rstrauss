#' Fit Strauss process to data using Logistic likelihood method
#' 
#' 
#' @param x a list with elements $x:data.frame of locations (ncol 2 or 3), $bbox: matrix with columns giving window bounds
#' @param R the fixed known interaction range. Use profile pl to infer this.
#' @param rho density of Poisson dummies
#' @export

fstrauss.logistic <- function(x, R, rho, ...) {
  bbox <- if(is.list(x)) x$bbox else apply(x, 2, range)
  x <- ( if(is.list(x)) x$x else x )
  N <- nrow(x)
  dim <- ncol(bbox)
  V <- prod(apply(bbox, 2, diff))
  #' dummy intensity
  if(missing(rho)) rho <- 4*N/V
  #'
  K <- rpois(1, rho*V)
  dummies <- apply(bbox, 2, function(ra) runif(K, ra[1], ra[2]))
  colnames(dummies) <- colnames(x)
  #'
  #' compute statistics
  loc <- rbind(x, dummies)
  nlist <- geom(loc, to=1:N, r=R)
  z <- rep(1:0, c(N, K))
  #'
  degs <- sapply(nlist, length)
  X <- cbind(degs) # intercept automatically
  #'
  bdry <- bbox_distance(loc, bbox) > R
  #' offset
  H <- rep(log(1/rho), K+N)
  #' fit
  fit <- glm(z~X, family="binomial", offset=H, subset=bdry)
  #'
  #'that's it
  coef  <- c(exp(fit$coef), R)
  names(coef) <- c("beta", "gamma", "r_given")
  list(theta=coef, fit=fit, logLik=logLik(fit))
}


