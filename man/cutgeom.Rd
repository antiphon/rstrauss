\name{cutgeom}
\alias{cutgeom}
\title{Cut the edges of from a geometric graph}
\usage{
cutgeom(x, nlist, r)
}
\arguments{
  \item{x}{point locations, n times d - matrix}

  \item{nlist}{neighbourhood list, n-long list with each
  element i a vector of neighbours of point x[i,]}

  \item{r}{the range at which to cut edges longer than it.}
}
\description{
Basic idea: Given a set of point locations, and some graph
on them, cut edges from list-of-vectors longer than range
r.
}

