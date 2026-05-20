# #' Fit bivariate Strauss process to data using Logistic likelihood method
# #' 
# #' 
# #' @param x a list with elements $x:data.frame of locations (ncol 2 or 3), $mark: class labels (factor with 2 levels), $bbox: matrix with columns giving window bounds
# #' @param R the fixed known interaction ranges, 3.
# #' @param rho densities of Poisson dummies, 2.
# #' @export
# 
# fstraussbiv.logistic2 <- function(x, R, rho, ...) {
#   if(is.null(x$mark)) stop("$mark element missing from x.")
#   z <- if(is.factor(x$mark)) x$mark else factor(x$mark)
#   if(nlevels(z)!=2 ) stop("factor(x$mark) should have two levels.")
#   if(length(R)!=3) stop("R needs to be vector of length 3, in the order c(R1, R2, R12)")
#   N <- table(z)
#   n <- sum(N)
#   z <- 1*(z==levels(z)[1])
#   i1 <- which(z==1)
#   i2 <- which(z==0)
#   bbox <- if(is.list(x)) x$bbox else apply(x, 2, range)
#   dim <- ncol(bbox)
#   V <- prod(apply(bbox, 2, diff))
#   loc <- ( if(is.list(x)) x$x else x )
#   
#   #' dummy intensity
#   if(missing(rho)) rho <- 4*N/V
#   #'
#   K <- sapply(rho*V, rpois, n=1)
#   dummies <- apply(bbox, 2, function(ra) runif(sum(K), ra[1], ra[2]))
#   z <- c(z, rep(1:0, K))
#   colnames(dummies) <- colnames(loc)
#   id1 <- which(z==1)
#   id2 <- which(z==0)
#   #'
#   #' compute statistics
#   u <- rbind(loc, dummies)
#   X <- cbind(z, 1-z)
#   #' 
#   #' interactions:
#   Rm <- max(R); im <- match(Rm, R)
#   nlists <- list()
#   nlists[[im]] <- geom(u, to=1:n, r=Rm)
#   for(i in setdiff(1:3, im))nlists[[i]] <- cutgeom(u, nlists[[im]], r=R[i])
#   #' covariates
#   #' interaction of type A
#   w11 <- 
#   X <- cbind(X, w11)
#   
#   #'
#   bdry <- bbox_distance(u, bbox) > max(R)
#   #' offset
#   H <- rep(log(1/rho), N+K)
#   #' obs vector: 1 for data, 0 for dummy
#   e <- rep(1:0, c(sum(N), sum(K)) )
#   #' fit
#   fit <- glm(e ~-1 + X, family="binomial", offset=H, subset=bdry)
#   #'
#   #'that's it
#   coef  <- c(exp(fit$coef))
#   names(coef) <- c("beta1", "beta2", "gamma1", "gamma2", "gamma12")
#   list(theta=coef, fit=fit, logLik=logLik(fit))
# }
# 
# 
