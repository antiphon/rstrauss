#' distance from points to bbox border
#'
#' @param x coordinate matrix
#' @param bbox bounding box 
#' @export

bbox_distance <- function(x, bbox){
  dim <- ncol(bbox)
  di <- function(loc){
    d <- sapply(1:dim, function(i) min(loc[i]-bbox[1,i], abs(loc[i]-bbox[2,i])))
    min(d)
  }
  c(apply(x, 1, di))
}

#' dilate bbox
#'
#' @param bbox bounding box 
#' @param r size of dilation
#' @export 
bbox_dilate <- function(bbox, r){
  apply(bbox, 2, function(ab) ab + c(-1,1)*r  )
}

#' Erode a 3d box
#' 
#' @export
bbox_erode <- function(bbox, r=0.1) {
  bbox[1,] <- bbox[1,] + r
  bbox[2,] <- bbox[2,] - r
  if(any(bbox[1,]>bbox[2,])) warning("Eroded too much.")
  bbox
}
