#' Mean field fit of Strauss process using logistic likelihood approximation
#' 
#' Using Jaakkola&Jordan 1996 quadratic variational approach to the Baddeley et al 2014 
#' logistic likelihood approximation of the point pattern likelihood.
#' 
#' @export

fstrauss.vbll <- function(x, R, rho, eps=1e-4, maxiter=100, verb=FALSE, ...){
  #stop("VB-LL not yet ready.")
  #
  dbb <- apply(x$bbox, 2, diff)
  V <- prod(dbb)
  dim <- ncol(x$bbox)
  N <- nrow(x$x)
  # dummy intensity
  if(missing(rho)) rho <- 4*N/V
  # dummy points
  K <- rpois(1, rho*V)
  dummies <- apply(x$bbox, 2, function(ra) runif(K, ra[1], ra[2]))
  colnames(dummies) <- colnames(x)
  
  cat2 <- if(verb) cat else function(...) NULL
  
  # neighbourhoods, sufficient statistic in Papangelou
  G  <- geom(rbind(x$x, dummies), to=1:N, r=R)
  S <- sapply(G, length)
  #
  # Compile standard glm formal data, inlclude intercept (1. column)
  X <- cbind(1,S)
  y <- c(rep(1, N), rep(0, K))
  offset <- log( 1/rho ) # same for all, data and dummies alike (@TODO: inhomogeneous version)
  #
  fit <- vblogit(y=y, X=X, offset=offset, verb=verb, ...)
  coef  <- c(exp(unname(fit$m[,1])), R)
  names(coef) <- c("beta", "gamma", "r_given")
  list(theta=coef, fit=fit, logLik=fit$logLik)
}

