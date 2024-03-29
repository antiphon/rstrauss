% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_LL-lj.R
\name{flennardjones.logistic}
\alias{flennardjones.logistic}
\title{Fit Lennard-Jones process to data using Logistic likelihood method}
\usage{
flennardjones.logistic(x, sigma0, rho, border_r, ...)
}
\arguments{
\item{x}{a list with elements $x:data.frame of locations (ncol 2 or 3), $bbox: matrix with columns giving window bounds}

\item{sigma0}{Initial guess of sigma. Used for numerical stability scaling.}

\item{rho}{density of Poisson dummies}

\item{border_r}{Border correction range. If missing, will use sigma0.}

\item{...}{passed on to glm}
}
\description{
supports also 3D.
}
