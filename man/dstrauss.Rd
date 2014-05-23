\name{dstrauss}
\alias{dstrauss}
\title{density of a strauss process}
\usage{
dstrauss(x, beta, gamma, range, method = "OT")
}
\arguments{
  \item{x}{list with x: coordinate matrix (col-wise), bbox:
  bounding box (col-wise ranges)}

  \item{beta}{beta parameter}

  \item{gamma}{gamma parameter}

  \item{range}{interaction range}

  \item{method}{how to approximate.

  Currently only method="OT" for Ogata&Tanemura 1981 is
  available, using numerical summation.}
}
\description{
Now only for cuboidal window
}

