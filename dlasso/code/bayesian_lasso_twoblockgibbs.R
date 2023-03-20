library(mvtnorm)
library(invgamma)
#library(GeneralizedHyperbolic)
library(GIGrvg)
library(statmod)
library(SuppDists)

my_rmvt <- function(n, mu, sigma, nu) {
  return( t(t(rmvnorm(n, rep(0, length(mu)), sigma) * sqrt(nu / rchisq(n, nu))) + mu))
}


bayesian_lasso_twoblockgibbs <- function(vy, mX, lambda, nsamples=10000)
{
  MAXITER <- nsamples
  
  lambda2 <- lambda*lambda
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  yTy <- t(vy)%*%vy
  
  vb <- rep(1,p)
  va <- rep(1,p)
  sigma2 <- 1
  a_til <- 1
  b_til <- 1
  
  mBeta   <- matrix(0,MAXITER,p)
  vsigma2 <- matrix(0,MAXITER,1)
  mB <- matrix(0,MAXITER,p)
  
  for (ITER in 1:MAXITER) 
  {
    mQ_inv <- XTX + lambda2*diag(1/va)
    mQ <- solve(mQ_inv,tol=1.0E-99)
    mSigma <- (b_til/a_til)*mQ
    vmu <- mQ%*%XTy
    vmu <- as.vector(vmu)
    
    vbeta <- my_rmvt(1, vmu, mSigma, 2*a + n)
    vbeta <- as.vector(vbeta)
    
    a_til <- a + 0.5*n
    b_til <- b + 0.5*(yTy - sum(vmu*XTy))
    b_til <- as.vector(b_til)
    
    sigma2 <- rinvgamma(1, shape=a_til, rate=b_til)
    sigma <- sqrt(sigma2)
    
    vnu <- sigma/(lambda*abs(vbeta))  
    #vB <- lambda2*vbeta*vbeta/sigma2 
    for (j in 1:p) {
      #va[j] <- rgig(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[i] <- rgig(1, lambda=0.5, chi=1, psi=vB[i])
      va[j] <- 1/rinvgauss(1, mean=vnu[j])
      #va[j] <- 1/rinvGauss(1, nu=vnu[j], lambda=1)
    }
    #vb <- as.vector(vb)
    va <- as.vector(va)
    
    
    mBeta[ITER,] <- vbeta
    vsigma2[ITER] <- sigma2
    #mB[ITER,] <- vb
  }
  
  return(list(mBeta=mBeta, vsigma2=vsigma2,mB=mB))
}
