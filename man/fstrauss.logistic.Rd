\name{fstrauss.logistic}
\alias{fstrauss.logistic}
\title{Fit Strauss process to data using Logistic likelihood method}
\usage{
fstrauss.logistic(x, R, rho)
}
\arguments{
  \item{x}{a list with elements $x:data.frame of locations
  (ncol 2 or 3), $bbox: matrix with columns giving window
  bounds}

  \item{R}{the fixed known interaction range. Use profile
  pl to infer this.}

  \item{rho}{density of Poisson dummies}
}
\description{
Fit Strauss process to data using Logistic likelihood
method
}

