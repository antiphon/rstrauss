#' Mean field fit of Strauss process using pseudolikelihood approximation
#' 
#' Not yet implemented properly.
#' 
#' @export

fstrauss.vb <- function(x, range, nx=25, eps=1e-4, maxiter=100, verb=FALSE){
  stop("VB-PL not yet ready.")
  dbb <- apply(x$bbox, 2, diff)
  V <- prod(dbb)
  abb <- dbb/dbb[1]
  n <- nrow(x$x)
  dim <- ncol(x$bbox)
  
  cat2 <- if(verb) cat else function(...) NULL
  #' grid
  ng <- round(abb * nx)
  gridc <- sapply(1:dim, 
                function(i) 
                  seq(x$bbox[1,i]+range, x$bbox[2,i]-range, length=ng[i]) )
  grid <- as.matrix(expand.grid(as.data.frame(gridc)))
  m <- nrow(grid)
  #' cell vol
  v <- V/prod(ng)
  #' neighbours
  G  <- geom(x$x, r=range)
  S <- sum(sapply(G, length))/2
  Gu <- geom(rbind(grid, x$x), to=(m+1):(n+m), from=1:m, r=range)
  Su <- sapply(Gu[1:m], length)
  #'
  iter <- 0
  #'
  priors <- list(gamma=c(1, 1), beta=c(1, 0.001))
  #'
  gamma0 <- 1
  beta0 <- n/V
  #'
  gseq <- seq(0, 1, length=500)
  post <- list(beta=c(priors$beta[1]+n, 0), gamma=(gseq*0+1)/length(gseq))
  bhist <- NULL
  ghist <- NULL
  loop <- TRUE
  while(loop){
    bhist<-c(bhist, b0 <- post$beta[2])
    g0 <- post$gamma
    #' update beta
    s <- sapply(Su, function(deg) sum( post$gamma * (gseq^deg) ) ) 
    post$beta[2] <- priors$beta[2] + v*sum(s)
    #' update gamma
    mub <- post$beta[1]/post$beta[2]
    
    gf <- function(g) 
      g^(priors$gamma[1]+S-1) * (1-g)^(priors$gamma[2]-1) * exp(-mub*v*sum(g^Su))
    
    post$gamma <- sapply(gseq, gf)
    post$gamma <- post$gamma/sum(post$gamma)
    
    # loop ending
    dif <- max( sum(abs(g0-post$gamma)), abs(b0-post$beta[2]))
    conv <- eps < dif 
    loop <- ( (iter<-iter+1) < maxiter )  & conv
    cat2("\r", iter," [mu_beta=", mub,"]        ")
    ghist <- c(ghist,   sum(gseq * post$gamma ))
  }
  cat2("\n")
  post$gamma <- cbind(gamma=gseq, density=post$gamma)
  post$stats <- list(iter=iter, bhist=bhist, ghist=ghist)
  #' likelihood f(x)
  
  
  post
}

