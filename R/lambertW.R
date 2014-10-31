#' Lambert's W function order 0
#' 
#' shamelessly copied from spatstat.
#' 
#' @export

LambertW <- function (x) 
{
  if (require(gsl, quietly = TRUE)) 
    return(gsl::lambert_W0(x))
  result <- rep.int(NA_real_, length(x))
  yexpyminusx <- function(y,x){y*exp(y)-x}
  for (i in which(is.finite(x) & (x >= 0))) result[i] <- uniroot(yexpyminusx, 
                                                                 c(0, x[i]), x = x[i])$root
  return(result)
}
