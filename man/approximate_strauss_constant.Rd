\name{approximate_strauss_constant}
\alias{approximate_strauss_constant}
\title{Approximate Strauss constant (master function)}
\usage{
approximate_strauss_constant(beta, gamma, range, bbox, method = "P", ...)
}
\arguments{
  \item{beta}{beta}

  \item{gamma}{gamma}

  \item{range}{range}

  \item{bbox}{bounding box}

  \item{method}{OT, PS or P}

  \item{...}{further details for methods

  Conditioning on number of points can be done with OT and
  P by setting the named parameter n.

  Methods:

  OT : Direct approx., Ogata&Tanemura 1981.
  \code{\link{approximate_strauss_constant_OT}}

  P: Direct approx., Penttinen 1998.
  \code{\link{approximate_strauss_constant_penttinen}}

  PS: MC unbiased estimate, Path sampling Berthelsen&Moller
  2003. \code{\link{approximate_strauss_constant_PS}}}
}
\description{
Approximate the normalizing constant (Z) in the Strauss
density.
}

