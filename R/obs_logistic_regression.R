# #' Fit logistic regression model using VB approximation
# #'
# #' @import Matrix 
# #' @export
# 
# vb.logit<-function(y, X, offset, eps=1e-2, m0, S0, S0i, xi0, verb=FALSE, maxiter=100, ...) {
#   ### Logistic regression using JJ idea. Ormeron00 notation.
#   ## p(y, w, t) = p(y | w) p(w | t) p(t) 
#   ##
#   ## Y ~ Bern(logit(Xw + offset))
#   ## w  ~ N(m0, S0) iid
#   ##
#   ## "*0" are fixed priors.
#   ##
#   cat2 <- if(verb) cat else function(...) NULL
#   
#   ## Write 
#   X <- Matrix(X)
#   N <- length(y)
#   K <- ncol(X)
#   #'
#   #'
#   #' offset
#   if(missing('offset')) offset <- 0
#   if(length(offset)<N) offset <- rep(offset, N)[1:N]
#   #
#   ## Priors and initial estimates.
#   if(missing(S0))S0   <- Diagonal(1e5, n=K)
#   if(missing(S0i))S0i <- solve(S0)
#   if(missing(m0))m0   <- rep(0, K)
#   if(missing(xi0))xi0   <-rep(1, N)
#   est <- list(m=m0, S=S0, Si=S0i, xi=xi0)
#   #'
#   #' Constants:
#   OO <- offset%*%t(offset)
#   LE_CONST <- -0.5*t(m0)%*%S0i%*%m0 -0.5*determinant(S0)$mod - 0.5*sum(OO)
#   Sm0 <- S0i%*%m0
#   #
#   ## helper functions needed:
#   lambda <- function(x)  -tanh(x/2)/(4*x)
#   gamma <- function(x)  x/2 - log(1+exp(x)) + x*tanh(x/2)/4
#   ## 
#   update <- function(v){
#     #' Lambda matrix
#     L <- Diagonal( x = lambda(v$xi) )
#     #' update post covariance
#     Si <- S0i - 2 * t(X)%*%L%*%X
#     S <- solve(Si)
#     #' update post mean
#     m <- S%*%( t(X)%*%( (y-0.5) + 2*L%*%offset ) + Sm0  )
#     #' update variational parameters
#     xi2 <- diag(  X%*%(S+m%*%t(m))%*%t(X) + OO + 2*(X%*%m)%*%t(offset) )
#     xi <- sqrt(xi2)
#     list(m=m, S=S, Si=Si, xi=xi, L=L)
#   }
#   
#   ###
#   ## loop
#   le <- -Inf
#   loop <- TRUE
#   le_hist <- -Inf
#   while(loop){
#     old <- le
#     # update all
#     est <- update(est)
#     ## compute the log evidence
#     le <- 0.5*determinant(est$S)$mod  + sum( gamma(est$xi) ) + 
#       t(offset)%*%est$L%*%offset + 0.5*t(est$m)%*%est$Si%*%est$m + LE_CONST    
#     #' check convergence
#     d <- abs(old - le)[1]
#     loop <- d > eps
#     le_hist <- c(le_hist, le[1])
#     cat2("diff:", d, "             \r")
#   }
#   cat2("\n")
#   ## done.
#   est$logLik <- le[1]
#   est$logp_hist <- le_hist
#   ## return
#   est
# }
