#' Compute the marginal likelihood for Bayesian Strauss inference
#' 

marginal_likelihood <- function(stats, R, bbox, priors, mciter=1000){
  stop("Not yet implemented.")
  samp <- function(n) {
    cbind(beta=rgamma(n, priors$beta[1], priors$beta[2]), 
          gamma=rbeta(n, priors$gamma[1], priors$gamma[2]))
  }
  pf <- function(bg) (stats[1]*log(bg[1])+stats[2]*log(bg[2])-approximate_strauss_constant(bg[1],bg[2], R, bbox))
  grid <- samp(1000)
  p <- apply(grid, 1, pf)
  mean(p, na.rm=T)
}

