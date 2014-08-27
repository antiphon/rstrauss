\name{fgibbs.logistic}
\alias{fgibbs.logistic}
\title{Fit exp-Gibbs processes to data using Logistic likelihood method}
\usage{
fgibbs.logistic(x, par, model, rho, ...)
}
\arguments{
  \item{x}{a list with elements $x:data.frame of locations
  (ncol 2 or 3), $bbox: matrix with columns giving window
  bounds}

  \item{par}{fixed parameters.}

  \item{model}{which model to use}

  \item{rho}{density of Poisson dummies}
}
\description{
models: "strauss", "areainteraction"
}

