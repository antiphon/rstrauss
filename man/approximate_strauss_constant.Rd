\name{approximate_strauss_constant}
\alias{approximate_strauss_constant}
\title{Approximate Strauss constant (master function)}
\usage{
approximate_strauss_constant(beta, gamma, range, bbox, method = "R", ...)
}
\arguments{
  \item{beta}{beta}

  \item{gamma}{gamma}

  \item{range}{range}

  \item{bbox}{bounding box}

  \item{method}{OT, PS, P or R. Default: R}

  \item{...}{further details for methods

  Conditioning on number of points can be done with OT and
  P by setting the named parameter n.

  Methods:

  OT : Direct approx., Ogata&Tanemura 1981.
  \code{\link{approximate_strauss_constant_OT}}

  P: Direct approx., Penttinen 1998.
  \code{\link{approximate_strauss_constant_penttinen}}

  PS: MC unbiased estimate, Path sampling Berthelsen&Moller
  2003. \code{\link{approximate_strauss_constant_PS}}

  R: My own version for unconditional Z. Based on a normal
  approximation (not published). Works when beta*Area >> 0}
}
\description{
Approximate the normalizing constant (Z) in the Strauss
density.
}

