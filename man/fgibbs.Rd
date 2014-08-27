\name{fgibbs}
\alias{fgibbs}
\title{Fit stationary and isotropic exp-family Gibbs process in 2D and 3D}
\usage{
fgibbs(x, par, model, approx = "LL", ...)
}
\arguments{
  \item{x}{Point pattern, see details for format.}

  \item{par}{Fixed parameters}

  \item{model}{interaction model.}

  \item{approx}{Which approximation method to use. Default:
  "LL". See details.}

  \item{...}{Passed on to the method.

  Data can be given in several formats: \describe{
  \item{matrix}{Columnwise coordinates. Bounding box
  computed.} \item{ppp}{2D or 3D pattern from spatstat.}
  \item{list}{Should have elements $x columnwise matrix of
  coordinates and $bbox bounding box (col-wise ranges).} }

  Three approximation methods are available: \describe{
  \item{LL}{Logistic likelihood approximation (Baddeley et
  al 2014). See \link{fstrauss.logistic}} \item{D}{Direct
  optimization using approximation of the normalizing
  constant. See \link{fstrauss.direct}} \item{B}{Bayesian
  version of D. \link{fstrauss.bayes}} \item{VB}{Bayesian
  fit using logistic likelihood and variational methods
  \link{fstrauss.vbll}} }}
}
\description{
NOTE: This is an early version.
}

