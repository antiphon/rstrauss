#' Profile pseudolikelihood for Strauss
#' 
#' Compute the likelihood profile for vector of interaction ranges. (basically a for-loop)
#' 
#' @param R vector of R, interaction range, values to try out
#' @param ... parameters passed to \\code{link{fstrauss}}
#' @export

fstrauss.profile <- function(R, ..., verb=FALSE) {
  coefs <- lpl <- NULL
  cat2 <- if(verb) cat else function(...) NULL
  for(i in 1:length(R)) {
     res <- fstrauss(..., R=R[i])
     coefs <- rbind(coefs, res$theta)
     lpl <- c(lpl, res$logLik)
     cat2(".")
  }
  cat2("\n")
  o <- which.max(lpl)
  theta_opt <- coefs[o,]
  names(theta_opt)[3] <- "r_profiled"
  lpl_opt <- lpl[o]
  list(opt=list(theta=theta_opt, logLik=lpl_opt), R=R, coefs=coefs, logLik=unname(lpl))
}