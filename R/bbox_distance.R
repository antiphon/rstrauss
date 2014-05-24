#' distance from points to bbox border
#'
#' @param x coordinate matrix
#' @param bbox bounding box 
#' @export
bbox_distance <- function(x, bbox){
  di <- function(loc){
    d <- sapply(1:ncol(bbox), function(i) min(loc[i]-bbox[1,i], abs(loc[i]-bbox[2,i])))
    min(d)
  }
  c(apply(x, 1, di))
}