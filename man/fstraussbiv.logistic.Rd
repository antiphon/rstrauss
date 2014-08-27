\name{fstraussbiv.logistic}
\alias{fstraussbiv.logistic}
\title{Fit bivariate Strauss process to data using Logistic likelihood method}
\usage{
fstraussbiv.logistic(x, R, rho, ...)
}
\arguments{
  \item{x}{a list with elements $x:data.frame of locations
  (ncol 2 or 3), $mark: class labels (factor with 2
  levels), $bbox: matrix with columns giving window bounds}

  \item{R}{the fixed known interaction ranges, 3.}

  \item{rho}{densities of Poisson dummies, 2.}
}
\description{
Fit bivariate Strauss process to data using Logistic
likelihood method
}

