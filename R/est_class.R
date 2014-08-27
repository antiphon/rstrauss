#' Fit result object
#' 
#' @export

print.straussfit <- function(x, ...) {
  cat("Method:", x$method, "\n")
  cat("Range:", x$theta[3], "\n")
  cat("Estimates:\n")
  print(x$theta[1:2])
}

#' Fit result object
#' 
#' @export

print.gibbsfit <- function(x, ...) {
  cat("Model:", x$model$model, "\n")
  cat("Given parameters:", paste(x$model$par, collapse=", "), "\n")
  cat("Method:", x$method, "\n")
  cat("Estimates:\n")
  print(x$theta)
}