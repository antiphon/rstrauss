#' Computes the sufficient statistics for exp-Gibbs models
#' 
#' @export

gibbs_pot <- function(x, par, model, ...) {
  models <- c("strauss", "lennardjones")
  i <- match(model, models)
  if(is.na(i)) stop(paste("'model' should be one of:", paste(models, collapse=", ")))
  pot <- get(paste0("gibbs_pot.", models[i]))
  pot(x, par, ...)
  
}

#' strauss
#' 
gibbs_pot.strauss <- function(x, par, ...) {
  nlist <- geom(x, r=par[1], ...)
  degs <- sapply(nlist, length)
  coef <- function(theta){
    theta <- as.data.frame(matrix(c(exp(theta), par[1]), nrow=1))
    names(theta) <- c("beta", "gamma", "range")
    theta
  }
  list(X=cbind(degs), border=par[1], coef=coef)
}

#' lennard-jones
#' 
gibbs_pot.lennardjones <- function(x, par, to) {
  d <- as.matrix(dist(x))
  diag(d) <- Inf
  d <- d[,to]
  #' rescale
  s0 <- 1
  s0 <- min(d[upper.tri(d, F)])
  if(s0==0) warning("minimum pairwise distance is 0.")
  d <- d/s0
  #' first
  v1 <- -apply(1/d^12, 1, sum)
  v2 <- apply(1/d^6, 1, sum)
  #' coefficient transformation
  cf <- function(theta) {
    theta[2] <- -theta[2] 
    eps <- theta[3]^2/(4*theta[2])
    sigma <- s0*(theta[3]/eps)^(1/6)
    list(canonical=data.frame(log_beta=theta[1], theta1=theta[2], theta2=theta[3], row.names=""),
         original=data.frame(beta=exp(theta[1]), eps=eps, sigma=sigma, row.names=""))
  }
  list(X=cbind(v1, v2), border=2.5*s0, coef=cf)
}




