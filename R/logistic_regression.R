#' Fit logistic regression model using VB
#'
#' @import Matrix 
#' @export

vb.logit<-function(y, X, offset, eps=1e-2, m0, S0, S0i, l0, verb=FALSE, maxiter=100,...) {
  ### Logistic regression using JaakkolaJordan96 approximation.
  ## p(y, w, t) = p(y | w) p(w | t) p(t) 
  ##
  ## y ~ Bern(logit(w^T x+offset))
  ## w  ~ N(0, t^-1) iid
  ## t  ~ Gamma(a0, b0)
  ##
  ## "*0" are fixed.
  ##
  ## Variational approximation q(w,t)=q(w)q(t) 
  ## And logistic function is approximated by
  ## logit(x) ~= logit(z)exp((x-z)/2 - lambda(z)(x^2-z^2))
  ## where lambda(z)=(logit(z)-0.5)/(2*z)
  ## and z is the variational parameter (auxiliary parameter) 
  ## from Taylor approximation.
  ##
  ##
  ## Write 
  cat2 <- if(verb) cat else function(...) NULL
  X <- Matrix(X)
  N <- length(y)
  K <- nrow(X)
  #y[y<1] <- -1
  D <- X%*%(y-0.5)
  #'
  #' offset
  if(missing('offset')) offset <- 0
  if(length(offset<N)) offset <- rep(offset, N)[1:N]
  #
  ## helper functions needed:
  logit <- function(x)1/(1+exp(-x))
  lambda <- function(l)  (logit(l)-0.5)/(2*l)
  CF <- function(x) sum( x/2 - log(1+exp(x)) + x*tanh(x/2)/4 )
  #
  #
  OO <- offset*offset
  SO <- 0#sum((y-0.5)*offset)
  ## Update formulas are no longer so
  ## easy computations:
  ## We udpate the w parameters (mN, tN) together with the variational parameter
  qw<-function(w) {
    le <- c(lambda(w$l))
    L <- Diagonal(n=N, x=le)
    Si <- S0i + 2*X%*%L%*%t(X)
    S <- solve(Si)
    off <- -2 * X%*%L%*%offset
    m <- S%*%( D + S0i%*%m0 + off)
    dim(m)
    V <- 2 * offset * t(X)%*%m
    l <- diag(t(X)%*%( S + m%*%t(m))%*%(X)) + V + OO
    l <- as.numeric(sqrt(l))
    list(m=m, S=S, Si=Si, l=l)
  }
  ## 
  ## initialize.
  if(missing(S0))S0   <- Diagonal(1000000, n=K)
  if(missing(S0i))S0i <- solve(S0)
  if(missing(m0))m0   <- rep(0, K)
  if(missing(l0))l0   <-rep(2,N)
  wN <- list(m=m0, S=S0, Si=S0i, l=l0)
  ## loop
  le<--Inf
  loop <- TRUE
  le_hist <- -Inf
  iter <- 0
  while(loop){
    old <- le
    wN <- qw(wN)
    lN <- wN$l
    ## compute the log evidence
    le <- as.numeric(0.5*determinant(wN$S)$modulus + 
                       0.5*t(wN$m)%*%wN$Si%*%wN$m + CF(wN$l) + SO)
    
    d <- abs(old - le)
    loop <- d > eps & ((iter <- iter+1) < maxiter)
    le_hist <- c(le_hist, le)
    cat2("d:", d, "             \r")
  }
  cat2("\n")
  ## done.
  if(iter == maxiter) warning("Didn't converge.")
  ## log evidence constant
  le_const <- as.numeric(- 0.5*determinant(S0)$modulus - 0.5*t(m0)%*%S0i%*%m0 )
  wN$logp <- as.numeric(le_hist) + le_const
  ## return
  wN
}
