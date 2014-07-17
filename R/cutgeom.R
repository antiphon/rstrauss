#' Cut the edges of from a geometric graph
#' 
#' Basic idea: Given a set of point locations, and some graph on them, cut edges from list-of-vectors longer than range r.
#' 
#' @param x point locations, n times d - matrix
#' @param nlist neighbourhood list, n-long list with each element i a vector of neighbours of point x[i,]
#' @param r the range at which to cut edges longer than it.
#' 
#' @export

cutgeom <- function(x, nlist, r) {
  x <- as.matrix(x)
  n <- nrow(x)
  if(length(nlist)!=n) stop("nlist does not match x.")
  new_nlist <- c_cutgeom(x, nlist, r)
  new_nlist
}