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
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  
  vb <- rep(1,p)
  va <- rep(1,p)
  sigma2 <- 10
  
  mBeta    <- matrix(0,MAXITER,p)
  vsigma2  <- matrix(0,MAXITER,1)
  vlambda2 <- matrix(0,MAXITER,1)
  mB <- matrix(0,MAXITER,p)
  mA <- matrix(0,MAXITER,p)
  
  for (ITER in 1:MAXITER) 
  {
    mQ_inv <- XTX + lambda2*diag(1/va)
    mQ <- solve(mQ_inv,tol=1.0E-99)
    mSigma <- sigma2*mQ
    vmu <- mQ%*%XTy
    vmu <- as.vector(vmu)
    vbeta <- t(rmvnorm(n = 1, mean = vmu, sigma = mSigma))
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
    }
    
    
    mBeta[ITER,] <- vbeta
    vsigma2[ITER] <- sigma2
    vlambda2[ITER] <- lambda2
    #mB[ITER,] <- vb
    mA[ITER,] <- va
    
    if ((ITER%%1000)==0)
      cat("ITER=",ITER,"\n")
  }
  
  return(list(mBeta=mBeta, vsigma2=vsigma2,vlambda2=vlambda2,mB=mB, mA=mA))
}
