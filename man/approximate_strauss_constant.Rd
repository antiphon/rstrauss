\name{approximate_strauss_constant}
\alias{approximate_strauss_constant}
\title{Master: Approximate Strauss constant}
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

  Methods: OT : Direct approx., Ogata&Tanemura 1981.
  \code{approximate_strauss_constant_OT}

  P: Direct approx., Penttinen
  1998\code{approximate_strauss_constant_penttinen}

  PS: MC unbiased estimate, Path sampling Berthelsen&Moller
  2003. \code{approximate_strauss_constant_PS}}
}
\description{
Master: Approximate Strauss constant
}

