# #' Mean field fit of Strauss process using logistic likelihood approximation
# #' 
# #' Using Jaakkola&Jordan 1996 quadratic variational approach to the Baddeley et al 2014 
# #' logistic likelihood approximation of the point pattern likelihood.
# #' 
# #' @export
# 
# fstraussbiv.vbll <- function(x, R, rho, eps=1e-4, maxiter=100, verb=FALSE, ...){
#   if(is.null(x$mark)) stop("$mark element missing from x.")
#   z <- if(is.factor(x$mark)) x$mark else factor(x$mark)
#   if(nlevels(z)!=2 ) stop("factor(x$mark) should have two levels.")
#   if(length(R)!=3) stop("R needs to be vector of length 3, in the order c(R1, R2, R12)")
#   N <- table(z)
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
#   #' to type 1
#   nlist11 <- geom(u, from=id1, to=i1, r=R[1])
#   X <- cbind(X, sapply(nlist11, length))
#   #' to type 2
#   nlist22 <- geom(u, from=id2, to=i2, r=R[2])
#   X <- cbind(X, sapply(nlist22, length))
#   #' 1-2
#   nlist12 <- geom(u, from=id1, to=i2, r=R[3])
#   nlist21 <- geom(u, from=id2, to=i1, r=R[3])
#   X <- cbind(X, sapply(nlist21, length)+sapply(nlist12, length))
#   #'
#   bdry <- bbox_distance(u, bbox) > max(R)
#   #' offset
#   H <- rep(log(1/rho), N+K)
#   #' obs vector: 1 for data, 0 for dummy
#   e <- rep(1:0, c(sum(N), sum(K)) )
#   #' fit
#   fit <- vb.logit(y=e, X=t(X), offset=H, verb=verb, ...)
#   coef  <- exp(unname(fit$m[,1]))
#   names(coef) <- c("beta1", "beta2", "gamma1", "gamma2", "gamma12")
#   list(theta=coef, fit=fit, logLik=fit$logp[length(fit$logp)])
# }
# 
