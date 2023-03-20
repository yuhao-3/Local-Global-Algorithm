calc_standard_errors <- function(vy, mX, lambda, vbeta, sigma2, vd, type=c("original","trimmed")) 
{
  lambda2 <- lambda*lambda
  sigma4 <- sigma2*sigma2 
  sigma6 <- sigma2*sigma4 
  
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  
  mH <- matrix(0,p+1,p+1)
  
  vy_hat <- mX%*%vbeta
  sigma2_hat <- sum(vy*(vy - vy_hat))/n
  
  mH[1:p,1:p] <-  -(XTX)/sigma2 
  mH[1:p,1+p] <-  -0.5*lambda2*(vd*vbeta)/sigma4
  mH[1+p,1:p] <- mH[1:p,1+p]
  mH[1+p,1+p] <- -(b + 0.5*n*sigma2_hat)/sigma6 + 0.25*lambda2*sum(vbeta*vd*vbeta)/sigma6
  
  if (type=="original") {
    mSigma <- solve(-mH,tol=1.0E-99)
    vse <- sqrt(diag(mSigma))
  } else {
    inds <- which(vbeta==0)
    mSigma <- solve(-mH[inds,inds],tol=1.0E-99)
    vse <- rep(0,p)
    vse[inds] <- sqrt(diag(mSigma))
  }
  
  return(list(mSigma=mSigma,vse=vse))
}







bayesian_lasso_em <- function(vy, mX, lambda, vbeta_init, sigma2_init, VERBOSE=FALSE, CALC_SE=TRUE)
{
  MAXITER <- 5000
  TOL <- 1.0E-6
  
  lambda2 <- lambda*lambda
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  yTy <- drop(t(vy)%*%vy)
  
  vbeta  <- vbeta_init
  sigma2 <- sigma2_init 
  vb <- rep(1,p)
  
  vtheta <- c(vbeta,sigma2)
  for (ITER in 1:MAXITER) 
  {
    # M-step
    vbeta <- solve(XTX + lambda2*diag(vb),tol=1.0E-99)%*%XTy
    
    #inds <- which(vb<1.0E8)
    #vbeta<- solve(XTX[inds,inds] + lambda2*diag(vb[inds]),tol=1.0E-99)%*%(XTy[inds])
    #vbeta[-inds] <- 0
    
    sigma2_hat <- (yTy - sum(XTy*vbeta))/n
    sigma2 <- (b + 0.5*n*sigma2_hat)/(a + 0.5*(n + p) + 1)
    sigma <- sqrt(sigma2)
    
    # E-step
    #vb[inds] <- sigma/(lambda*abs(vbeta[inds]))
    #vb[-inds] <- 1.0E8
    #vb <- as.vector(vb)
    
    EPS <- 1.0E-14
    vb <- sigma/(lambda*sqrt(vbeta^2 + EPS))
    vb <- as.vector(vb)
    
    vtheta_old <- vtheta
    vtheta <- c(vbeta,sigma2)
    err <- max(abs(vtheta - vtheta_old))
    if (VERBOSE) {
      cat("ITER=",ITER,"err=",err,"sigma2=",sigma2,"\n") #vbeta=",vbeta,"\n")
    }
    if (err < TOL) {
      break;
    }
  }
  
  # Set very small estimates to zero
  inds <- which(abs(vbeta)<1.0E-5)
  vbeta[inds] <- 0
  
  if (CALC_SE) {
    res <- calc_standard_errors(vy, mX, lambda, vbeta, sigma2, vb, type="original")
  }
  else {
    res <-list()
    res$mSigma <- NULL
    res$vse <- NULL
  }
  
  return(list(vbeta=vbeta,sigma2=sigma2, vb=vb, mSigma=res$mSigma, vse=res$vse))
}

