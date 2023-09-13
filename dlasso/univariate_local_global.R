
source("lasso_distribution.R")


univariate_local_global <- function(vy, mX, lambda, vmu_init, mSigma_init, a_til, b_til)
{
  ## Initialization
  MAXITER = 500
  TOL = 1.0E-7
  n = nrow(mX)
  p = ncol(mX)
  
  XTX = t(mX)%*%mX
  XTy = t(mX)%*%vy
  yTy = t(vy)%*%vy
 
  ## Initial global mean and covariance
  vmu_glob = vmu_init
  if(p<n){
    mSigma_glob = mSigma_init
  }else{
    mSigma_glob = diag(mSigma_init)
  }
  
  
  # Later will be changed
  c_val = lambda * (exp(lgamma(a_til+0.5) - lgamma(a_til) - 0.5*log(b_til)))
  
  E_invsigma2 <- a_til/b_til 
  
  # Record the local parameter
  va = c()
  vb = c()
  vZ = c()
  
  # Initial value Local parameters
  vtheta <- c(vmu_glob)
  
  for (ITER in 1:MAXITER)
  {
    ## Local Update
    for (j in 1:p)
    {
      ##########################################################################
      
      # Local Update
      
      # Store parameter for previous iteration
      vmu_old = vmu_glob
      mSigma_old = mSigma_glob
      
      vt = matrix(mSigma_glob[-j,j])/mSigma_glob[j,j]
      vs = matrix(vmu_glob[-j]) - vt*vmu_glob[j]
      
      ## Update local parameters
      a_val = E_invsigma2*(XTX[j,j] + XTX[j,-j]%*%vt)
      b_val = E_invsigma2*(XTy[j] - XTX[j,-j]%*%vs)
      
      
      ## Calculate Local mean and variance
      mu_star = elasso(a_val,b_val,c_val)
      sigma2_star = vlasso(a_val,b_val,c_val)
      
      mu_star <- as.vector(mu_star)
      sigma2_star <- as.vector(sigma2_star)
      
      # Record local parameter
      va[j] = a_val
      vb[j] = b_val
      vZ[j] = zlasso(a_val,b_val,c_val)
      
      ##########################################################################
      
      # Global Update
      ## Update Mean
      
      vmu_glob[j] = mu_star
      vmu_glob[-j] = vmu_glob[-j] + mSigma_glob[j,-j]*(mu_star - vmu_old[j])/mSigma_glob[j,j]
      
      ## Update Covariance
      mSigma_glob[j,j] = sigma2_star
      mSigma_glob[j,-j] = (sigma2_star/mSigma_old[j,j])*matrix(mSigma_old[j,-j])
      mSigma_glob[-j,j] = t(mSigma_glob[j,-j])
      mSigma_glob[-j,-j] = matrix(mSigma_old[-j,-j],p-1,p-1) + matrix(mSigma_old[-j,j])%*%((sigma2_star - matrix(mSigma_old[j,j]))/(mSigma_old[j,j]^2))%*%t(matrix(mSigma_old[j,-j]))
            
      # Damping
      # mSigma_adjust = rho* mSigma_adjust + (1-rho)*mSigma_old
    }
    
    # Update Theta
    vtheta_old <- vtheta
    vtheta <- c(vmu_glob)
    
    ## Check stopping criterion    
    err <- max(abs(vtheta - vtheta_old))
    #cat("ITER=",ITER,"err=",err,"\n")
    if (err < TOL) {
      break;
    }
  }
  # Return Local and Global parameter
  return(list("vmu_til" = vmu_glob,"mSigma_til" = mSigma_glob, "a" = va, "b" = vb,"c_val"= c_val ,"Z" = vZ))
  
}
