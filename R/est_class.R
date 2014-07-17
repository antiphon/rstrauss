#' Fit result object
#' 
#' @export

print.straussfit <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Range:", x$theta[3], "\n")
  cat("Estimates:\n")
  print(x$theta[1:2])
}