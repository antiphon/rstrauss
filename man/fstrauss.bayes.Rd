\name{fstrauss.bayes}
\alias{fstrauss.bayes}
\title{Bayesian inference for stationary Strauss process}
\usage{
fstrauss.bayes(x, R, priors = list(beta = c(1, 1e-05), gamma = c(1, 1)),
  lower = c(1e-09, 1e-09), upper = c(Inf, 1), init, ...)
}
\arguments{
  \item{x}{a list with elements $x:data.frame of locations
  (ncol 2 or 3), $bbox: matrix with columns giving window
  bounds}

  \item{R}{the fixed known interaction range. Use profile
  pl to infer this.}

  \item{priors}{A list of priors. See details.}

  \item{lower}{Lower limits for searching (beta,gamma).}

  \item{upper}{Upper limits for searching (beta,gamma).}

  \item{init}{Initial values for optim.}

  \item{...}{passed on to
  \code{\link{approximate_strauss_constant}}

  Reduced sample border correction with radius R.

  We use Gamma(a,b) prior for beta and Beta(c,d) prior for
  Gamma (funny...) so priors-parameter should be a list
  with elements $gamma=c(a,b) and $beta=c(c,d).

  At the moment only the MAP is computed and returned in
  the $par element.

  See \link{approximate_strauss_constant} for approximation
  options.}
}
\description{
Using direct approximation of the normalizing constant.
}

