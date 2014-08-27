#' Interpret data
#' 
#' @param x data to be converted.
#' @export

convert_to_pp <- function(x){
  cls <- is(x)
  if("matrix"%in%cls | "data.frame"%in%cls){
    l <- x
    bbox <- apply(cls, 2, range)
  }
  else if("list"%in% cls){
    l <- x$x
    bbox <- if(is.null(x$bbox)) apply(cls, 2, range) else x$bbox
  }
  else if("ppp"%in%cls){
    l <- cbind(x$x, x$y)
    bbox <- cbind(x$window$xrange, x$window$yrange) 
  }
  else{
    stop("Can't interpret data x.")
  }
  list(x=l, bbox=bbox, dim=ncol(bbox), n=nrow(l))
}
