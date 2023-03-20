bayesian_lasso_mfvb <- function(vy, mX, lambda, sigma2_hat) 
{
  MAXITER <- 500
  TOL <- 1.0E-8
  
  lambda2 <- lambda*lambda
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  
  vmu_til  <- matrix(0,p,1)
  mSigma_til <- diag(1,p)
  
  a_til <- a + 0.5*(n + p)
  b_til <- sigma2_hat*a_til
  
  vd <- rep(1,p)
  
  vtheta <- c(vmu_til,mSigma_til, a_til, b_til)
  for (ITER in 1:MAXITER) 
  {
    # M-step
    mQ_inv <- XTX + lambda2*diag(vd)
    mQ <- solve(mQ_inv,tol=1.0E-99)
    vmu_til    <- mQ%*%XTy
    mSigma_til <- (b_til/a_til)*mQ  
    vsigma2_til <- diag(mSigma_til)
    
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vmu_til)^2) + 0.5*lambda2*sum(vd*vmu_til*vmu_til) + 0.5*sum(diag(mQ_inv%*%mSigma_til))
    
    # E-step
    vd <- sqrt(b_til/(a_til*(vmu_til^2 + vsigma2_til)))/lambda
    vd <- as.vector(vd)
    
    vtheta_old <- vtheta
    vtheta <- c(vmu_til,mSigma_til,a_til, b_til)
    err <- max(abs(vtheta - vtheta_old))
    cat("ITER=",ITER,"err=",err,"\n")  
    if (err < TOL) {
      break;
    }
  }
  
  
  
  return(list(vmu_til=vmu_til, mSigma_til=mSigma_til, a_til=a_til, b_til=b_til, vd=vd))
}