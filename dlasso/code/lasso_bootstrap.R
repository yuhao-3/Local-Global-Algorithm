
lasso_naive_bootstrap <- function(vy,mX,lambda,B=5000,VERBOSE=FALSE)
{
  n <- nrow(mX)
  p <- ncol(mX)
  mBeta <- matrix(0,B,p)
  for (i in 1:B) 
  {
    inds_boot <- sample(1:n,size=n,replace=TRUE)
    vy_boot <- vy[inds_boot]
    mX_boot <- mX[inds_boot,]
    
    res_boot <-
      bayesian_lasso_em(
        vy_boot,
        mX_boot,
        lambda,
        vbeta_init = rep(0,p),
        sigma2_init = 1,
        VERBOSE = FALSE,
        CALC_SE = FALSE
      )
    
    if (VERBOSE) {
      if ((i%%1000)==1) {
        cat("i=",i,"\n")
      }
    }
    
    mBeta[i,] <- res_boot$vbeta
  }
  
  return(mBeta)
}


lasso_modified_bootstrap <- function(vy,mX,lambda,B=5000,a=0.5,VERBOSE=FALSE)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  res <-
    bayesian_lasso_em(
      vy,
      mX,
      lambda,
      vbeta_init = rep(0,p),
      sigma2_init = 1,
      VERBOSE = FALSE,
      CALC_SE = FALSE)
  
  # Lasso estimator
  vbeta_hat <- res$vbeta
  
  # Modified estimator
  vbeta_til <- vbeta_hat*(abs(vbeta_hat)>a)
  
  # Modified residuals
  vr <- vy - mX%*%vbeta_til

  # Centered modified residuals
  ve <- vr - mean(vr)

  mBeta <- matrix(0,B,p)
  for (i in 1:B) 
  {
    inds_boot <- sample(1:n,size=n,replace=TRUE)
    
    ve_star <- ve[inds_boot]
    vy_star <- mX%*%vbeta_til + ve_star 
      
    res_boot <-
      bayesian_lasso_em(
        vy_star,
        mX,
        lambda,
        vbeta_init = rep(0,p),
        sigma2_init = 1,
        VERBOSE = FALSE,
        CALC_SE = FALSE)
    
    if (VERBOSE) {
      if ((i%%1000)==1) {
        cat("i=",i,"\n")
      }
    }
    
    mBeta[i,] <- res_boot$vbeta
  }
  
  return(mBeta)
}


