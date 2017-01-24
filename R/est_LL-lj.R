#' Fit Lennard-Jones process to data using Logistic likelihood method
#' 
#' supports also 3D.
#' 
#' @param x a list with elements $x:data.frame of locations (ncol 2 or 3), $bbox: matrix with columns giving window bounds
#' @param sigma0 Initial guess of sigma. Used for numerical stability scaling.
#' @param rho density of Poisson dummies
#' @param border_r Border correction range. If missing, will use sigma0.
#' @param ... passed on to glm
#' @export

flennardjones.logistic <- function(x, sigma0, rho, border_r, ...) {
  if(missing(sigma0)) stop("sigma0 missing.")
  if(missing(border_r)) border_r <- sigma0
  bbox <- if(is.list(x)) x$bbox else apply(x, 2, range)
  x <- ( if(is.list(x)) x$x else x )
  N <- nrow(x)
  dim <- ncol(bbox)
  V <- bbox_volume(bbox)
  # dummy intensity
  if(missing(rho)) rho <- 4*N/V
  #
  K <- rpois(1, rho*V)
  dummies <- apply(bbox, 2, function(ra) runif(K, ra[1], ra[2]))
  colnames(dummies) <- colnames(x)
  #
  # compute statistics
  
  loc <- rbind(x, dummies)
  D <- as.matrix(dist(loc)) / sigma0
  diag(D) <- Inf
  DS <- D[,1:N]
  #
  D6<-(1/DS)^6
  D12 <- D6^2
  # stability, like in spatstat
  D12[DS < .25] <- D6[DS < .25] <- -Inf
  D12[DS > 4]   <- D6[DS > 4] <- 0
  t2 <-   rowSums(D6)
  t1 <-   -rowSums(D12)
  #
  z <- rep(1:0, c(N, K))
  #
  X <- cbind(t1, t2) # intercept automatically
  #
  bdry <- bbox_distance(loc, bbox) > border_r & is.finite(t1) & is.finite(t2)
  # offset
  H <- rep(log(1/rho), K+N)
  # fit
  fit <- glm(z~X, family="binomial", offset=H, subset=bdry, ...)
  #
  # that's it
  # Backtransform
  co  <- fit$coef
  beta <- exp(co[1])
  theta <- c(beta, (co[2]/co[3])^(1/6)*sigma0,
             co[3]^2/(4*co[2]))
  names(theta) <- c("beta", "sigma", "epsilon")
  list(theta=theta, fit=fit, logLik=logLik(fit))
}


