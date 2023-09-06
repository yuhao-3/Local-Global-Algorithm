library(mvtnorm)
library(invgamma)
#library(GeneralizedHyperbolic)
library(GIGrvg)
library(statmod)
library(SuppDists)

bayesian_lasso_gibbs <- function(vy, mX, lambda, nsamples=10000)
{
  MAXITER <- nsamples
  
  
  if (is.numeric(lambda)) {
    lambda2 <- lambda*lambda
    USEPRIOR <- FALSE
  } else {
    if (is.list(lambda)) {
      USEPRIOR <- TRUE
      lambda_list <- lambda 
      lambda <- lambda_list$lambda_start
      lambda2 <- lambda*lambda
    }
  }
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  if (p<n) {
    XTX <- t(mX)%*%mX 
  } else {
    mI_n <- diag(rep(1, n))
  }
  XTy <- t(mX)%*%vy
  yTy <- t(vy)%*%vy
  
  vb <- rep(1,p)
  va <- rep(1,p)
  sigma2 <- 10
  
  mBeta    <- matrix(0,MAXITER,p)
  vsigma2  <- matrix(0,MAXITER,1)
  vlambda2 <- matrix(0,MAXITER,1)
  mB <- matrix(0,MAXITER,p)
  mA <- matrix(0,MAXITER,p)
  
  va_til <- matrix(0,MAXITER,1)
  vb_til <- matrix(0,MAXITER,1)
  vu_til <- matrix(0,MAXITER,1)
  vv_til <- matrix(0,MAXITER,1)
  mM <- matrix(0,MAXITER,p)
  mV <- matrix(0,MAXITER,p)
  
  for (ITER in 1:MAXITER) 
  {
    if (p<n) {
      mQ_inv <- XTX + lambda2*diag(1/va)
      mQ <- solve(mQ_inv,tol=1.0E-99)
      mSigma <- sigma2*mQ
      vmu <- mQ%*%XTy
      vmu <- as.vector(vmu)
      vbeta <- t(rmvnorm(n = 1, mean = vmu, sigma = mSigma))
    } else {
      vd  <- lambda2/va
      vd_inv <- 1/vd
      vu <- rnorm(p)*sqrt(vd_inv)
      vv <- mX%*%vu + rnorm(n)
      
      # Upper triangular Cholesky factorization
      mU <- chol(mX%*%diag(vd_inv)%*%t(mX) + mI_n)
      
      vw <- backsolve(mU, backsolve(t(mU), -vv, upper.tri=FALSE))
      vmu <- t(mX)%*%backsolve(mU, backsolve(t(mU), vy, upper.tri = FALSE))
      vmu <- vd*vmu
 
      vbeta <- vmu + sqrt(sigma2)*(vu + vd_inv*(t(mX)%*%vw))
    }
    vbeta <- as.vector(vbeta)
    
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vbeta)^2) + 0.5*lambda2*sum(vbeta*vbeta/va)
    b_til <- as.vector(b_til)
    
    sigma2 <- rinvgamma(1, shape=a_til, rate=b_til)
    #sigma <- sqrt(sigma2)
    
    #vnu <- sigma/(lambda*abs(vbeta))  
    vB <- lambda2*vbeta*vbeta/sigma2 
    for (j in 1:p) {
      va[j] <- rgig(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[i] <- rgig(1, lambda=0.5, chi=1, psi=vB[i])
      #vb[j] <- rinvgauss(1, mean=vmu[j])
      #vb[j] <- rinvGauss(1, nu=vnu[j], lambda=1)
    }
    #vb <- as.vector(vb)
    va <- as.vector(va)
    
    if (USEPRIOR) {
      if (lambda_list$prior=="gamma") {
        u <-   lambda_list$u
        v <-   lambda_list$v
        u_til <- u + 0.5*p
        v_til <- v + 0.5*sum(vbeta*vbeta/va)/sigma2
        lambda2 <- rgamma(1, shape=u_til, rate=v_til)
      } 
    } else {
      u_til <- NA
      v_til <- NA
    }
    
    
    mBeta[ITER,] <- vbeta
    vsigma2[ITER] <- sigma2
    vlambda2[ITER] <- lambda2
    #mB[ITER,] <- vb
    mA[ITER,] <- va
    
    va_til[ITER] <- a_til
    vb_til[ITER] <- b_til
    vu_til[ITER] <- u_til
    vv_til[ITER] <- v_til
    mM[ITER,] <- vmu
    
    if (p<n) {
      mV[ITER,] <- diag(mSigma)
    } else {
      mV[ITER,] <- vsigma_chol
    }
    
    if ((ITER%%1000)==0)
      cat("ITER=",ITER,"\n")
  }
  
  return(
    list(
      mBeta = mBeta,
      vsigma2 = vsigma2,
      vlambda2 = vlambda2,
      mB = mB,
      mA = mA,
      va_til = va_til,
      vb_til = vb_til,
      vu_til = vu_til,
      vv_til = vv_til,
      mM = mM,
      mV = mV
    )
  )
}



density_rw_sigma2 <- function(x, va_til, vb_til) 
{
  y <- x*0
  N <- length(va_til)
  for (i in 1:N) {
    y <- y + dinvgamma(x, shape=va_til[i], rate=vb_til[i])    
  }
  y <- y/N
  return(list(x=x,y=y))
}

density_rw_lambda2 <- function(x, vu_til, vv_til) 
{
  y <- x*0
  N <- length(vu_til)
  for (i in 1:N) {
    y <- y + dgamma(x, shape=vu_til[i], rate=vv_til[i])  
  }
  y <- y/N
  return(list(x=x,y=y))
}


density_rw_beta <- function(x, mM, mV,j) 
{
  y <- x*0
  N <- nrow(mM)
  for (i in 1:N) {
    y <- y + dnorm(x, mM[i,j], sqrt(mV[i,j]))  
  }
  y <- y/N
  return(list(x=x,y=y))
}





  
  
  