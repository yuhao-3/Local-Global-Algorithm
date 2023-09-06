
bayesian_lasso_mfvb <- function(vy, mX, lambda, sigma2_hat, a, b, maxiter, tol) 
{
  MAXITER <- maxiter
  TOL <- tol
  
  # Get model constants
  n <- nrow(mX)
  p <- ncol(mX)
  if (p<n) {
    XTX <- t(mX)%*%mX 
  } else {
    mI_n <- diag(rep(1, n))
  }
  XTy <- t(mX)%*%vy
  
  # Handle prior
  if (is.numeric(lambda)) {
    lambda2_til <- lambda*lambda
    USEPRIOR <- FALSE
  } else {
    if (is.list(lambda)) {
      USEPRIOR <- TRUE
      lambda_list <- lambda 
      u <-   lambda_list$u
      v <-   lambda_list$v
      lambda <- lambda_list$lambda_start
      lambda2_til <- lambda*lambda
      u_til <- u + 0.5*p
      v_til <- u_til/lambda2_til
    }
  }
  
  # Set variational parameters
  vmu_til  <- matrix(0,p,1)
  mSigma_til <- NULL
  a_til <- a + 0.5*(n + p)
  b_til <- sigma2_hat*a_til
  
  # Expected values of parameters under q
  sigma2inv_til <- 1/sigma2_hat
  va_til <- rep(1,p)

  # Set parameters to be monitored for convergence  
  vtheta <- c(vmu_til)
  for (ITER in 1:MAXITER) 
  {
    if (p<n) {
      # p<n case
      mQ_inv <- XTX + lambda2_til*diag(va_til)
      mQ <- solve(mQ_inv,tol=1.0E-99)
      vmu_til    <- mQ%*%XTy
      mSigma_til <- (1/sigma2inv_til)*mQ  
      vsigma2_til <- diag(mSigma_til)
    } else {
      # p>n case - via SMW
      vd <- lambda2_til*va_til
      vd_inv <- 1/vd
      
      XD <- t(t(mX)*vd_inv);
      #XD2 <- mX%*%diag(vd_inv)
      #print(max(abs(XD - XD2)))
      #ans <- readline()
      
      mQ_inv <- XD%*%t(mX) + mI_n
      mQ  <- solve(mQ_inv,tol=1.0E-99)
      QXD <- mQ%*%XD; 
      vmu_til <- t(QXD)%*%vy;
      vv <- t(XD*QXD)%*%matrix(1,n,1);
      vv <- as.vector(vv)
      vsigma2_til = (vd_inv - vv)/sigma2inv_til;
      
      #mU <- chol(mG)
      #mM3 <- backsolve(t(mU), mX, upper.tri=FALSE)
      #mM2 <- mM3*mM3
      #vv1 <- matrix(1,1,n)%*%mM2
      #vsigma2_til <- (1/Eq_sigma2inv)*(vd_inv - vd_inv*vd_inv*as.vector(vv1))
      #vmu_til <- t(mX)%*%backsolve(mU, backsolve(t(mU), vy, upper.tri = FALSE))
      #vmu_til <- vd_inv*vmu_til
    }
    
    #print(vmu_til)
    #ans <- readline()
    
    #print( cbind(vmu_til, vmu_til2, vsigma2_til, vsigma2_til2) )
    #ans <- readline()
    
    vmu2_til = vmu_til*vmu_til
    sum_Ea_vmu2 = sum(va_til*vmu2_til);
    
    b_til <- b + 0.5*(sum((vy - mX%*%vmu_til)^2) + lambda2_til*sum_Ea_vmu2 + p/sigma2inv_til) # simplified
    sigma2inv_til <- a_til/b_til

    
    if (USEPRIOR) {
      if (lambda_list$prior=="gamma") {
        v_til <- v + 0.5*sigma2inv_til*(sum_Ea_vmu2 + sum(va_til*vsigma2_til))
        lambda2_til <- u_til/v_til
      } 
    } else {
      v_til <- NA
    }
    
    va_til <- sqrt(1/(sigma2inv_til*lambda2_til*(vmu_til*vmu_til + vsigma2_til))) 
    va_til <- as.vector(va_til)
    
    vtheta_old <- vtheta
    vtheta <- c(vmu_til)
    err <- sqrt(sum(abs(vtheta - vtheta_old)^2))
    #cat("ITER=",ITER,"err=",err,"\n")  
    if (err < TOL) {
      break;
    }
  }
  
  return(
    list(
      vmu_til = vmu_til,
      mSigma_til = mSigma_til,
      vsigma2_til = vsigma2_til,
      a_til = a_til,
      b_til = b_til,
      u_til = u_til,
      v_til = v_til,
      sigma2inv_til = sigma2inv_til,
      lambda2_til = lambda2_til,
      va_til = va_til
    )
  )
}