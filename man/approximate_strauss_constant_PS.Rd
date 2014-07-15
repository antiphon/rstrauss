\name{approximate_strauss_constant_PS}
\alias{approximate_strauss_constant_PS}
\title{approximate Strauss constant using Path sampling}
\usage{
approximate_strauss_constant_PS(beta, gamma, range, bbox, nsim = 10,
  steps = 10, perfMaxIter = 1e+05, MHiter = 1000, ...)
}
\arguments{
  \item{nsim}{number of sims for MC expectation}

  \item{steps}{grid reso for gamma}

  \item{perfMaxIter}{max iter for the first perfect
  pattern}

  \item{MHiter}{iterations per MH simulations starting from
  perfect pattern}

  \item{...}{passed on to rstrauss simulator}
}
\value{
Value is the natural log of Z(theta)
}
\description{
approximate Strauss constant using Path sampling
}

