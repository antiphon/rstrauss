\name{approximate_strauss_constant_PS}
\alias{approximate_strauss_constant_PS}
\title{approximate Strauss constant using Path sampling}
\usage{
approximate_strauss_constant_PS(beta, gamma, range, bbox, nsim = 10,
  steps = 10, ...)
}
\arguments{
  \item{nsim}{number of sims for MC expectation}

  \item{steps}{grid reso for gamma}

  \item{...}{passed on to rstrauss simulator}
}
\value{
Value is the natural log of Z(theta)
}
\description{
approximate Strauss constant using Path sampling
}

