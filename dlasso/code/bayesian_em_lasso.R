################################################################################

calc_standard_errors <- function(vy, mX, lambda, vbeta, sigma2, vd, 
                                 type=c("original","trimmed")) 
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

################################################################################

bayesian_lasso_em <- function(vy, mX, lambda, a, b, 
                                   vbeta_init, sigma2_init, 
                                   trunc=0, VERBOSE=FALSE, CALC_SE=FALSE)
{
  MAXITER <- 5000
  TOL <- 1.0E-6
  BIG <- 1.0E6
  TINY1 <- 1.0E-4
  TINY2 <- 1.0E-12
  
  lambda2 <- lambda*lambda
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  if (p<n) {
    XTX <- t(mX)%*%mX 
  } else {
    mI_n <- diag(rep(1, n))
  }
  XTy <- t(mX)%*%vy
  yTy <- drop(t(vy)%*%vy)
  
  vbeta  <- vbeta_init
  sigma2 <- sigma2_init 
  sigma <- sqrt(sigma2)
  
  vb <- rep(1,p)
  #vb[inds] <- sigma/(lambda*abs(vbeta[inds]))
  #vb[-inds] <- BIG
  
  vtheta <- c(vbeta,sigma2)
  for (ITER in 1:MAXITER) 
  {
    inds <- which(abs(vbeta)!=0)
    s <- length(inds)
    
    if (s<n) {
      mX1 <- mX[,inds]
      XTX <- matrix(t(mX1)%*%mX1,s,s)
      if (length(inds)==1) {
        vbeta[inds] <- (XTX + lambda2*vb[inds])/XTy[inds]
      } else {
        vbeta[inds]  <- solve(XTX + lambda2*diag(vb[inds]),XTy[inds],tol=1.0E-99)
      }
      vbeta[-inds] <- 0
      
    } else {
      # M-step
      if (p<n) {
        vbeta <- solve(XTX + lambda2*diag(vb),XTy,tol=1.0E-99)
      } else {
        vd  <- lambda2*vb
        vd_inv <- 1/vd
        XD <- t(t(mX)*vd_inv);
        mQ_inv <- XD%*%t(mX) + mI_n
        mQ  <- solve(mQ_inv,tol=1.0E-99)
        QXD <- mQ%*%XD; 
        vbeta <- t(QXD)%*%vy;
      }
    }
    
    sigma2_hat <- (yTy - sum(XTy*vbeta))/n
    sigma2 <- (b + 0.5*n*sigma2_hat)/(a + 0.5*(n + p) + 1)
    sigma <- sqrt(sigma2)
    
    inds <- which(abs(vbeta)>TINY1)
 
    vbeta[-inds] <- 0
    vb[inds] <- sigma/(lambda*abs(vbeta[inds]+TINY2))
    vb[-inds] <- BIG
    vb <- as.vector(vb)
    
    vtheta_old <- vtheta
    vtheta <- c(vbeta,sigma2)
    err <- max(abs(vtheta - vtheta_old))
    if (VERBOSE!=0) {
      if ((ITER%%VERBOSE)==0) {
        cat("ITER=", ITER,"nnz(beta)=",s,"err=",err,"sigma2=",sigma2,"\n") #vbeta=",vbeta,"\n")
      }
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


################################################################################

bayesian_lasso_em_truncate <- function(vy, mX, lambda, 
                                       a, b, vbeta_init, sigma2_init, 
                                       VERBOSE=FALSE, CALC_SE=FALSE)
{
  MAXITER <- 5000
  TOL <- 1.0E-6
  BIG <- 1.0E6
  
  lambda2 <- lambda*lambda
  
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
    inds <- which(vb<BIG)
    vbeta[inds]  <- solve(XTX[inds,inds] + lambda2*diag(vb[inds]),XTy[inds],tol=1.0E-99)
    vbeta[-inds] <- 0
    
    sigma2_hat <- (yTy - sum(XTy*vbeta))/n
    sigma2 <- (b + 0.5*n*sigma2_hat)/(a + 0.5*(n + p) + 1)
    sigma <- sqrt(sigma2)
    
    # E-step
    vb[inds] <- sigma/(lambda*abs(vbeta[inds]))
    vb[-inds] <- BIG
    vb <- as.vector(vb)
    
    vtheta_old <- vtheta
    vtheta <- c(vbeta,sigma2)
    err <- max(abs(vtheta - vtheta_old))
    if (VERBOSE) {
      cat("ITER=",ITER, "nnz(beta)=",length(inds), "err=",err,"sigma2=",sigma2,"\n") #vbeta=",vbeta,"\n")
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
