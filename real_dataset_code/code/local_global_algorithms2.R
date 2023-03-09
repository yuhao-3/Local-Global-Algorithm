###############################################################################

## Local and Global Algorithm Two(Update each pair of variables at a time) 
## Expect an more accurate but slow performance
local_global_algorithm_2 <- function(vy, mX, lambda, params)
{
  ## Initialization
  MAXITER = 500
  TOL = 1.0E-8
  n = nrow(mX)
  p = ncol(mX)
  XTX = t(mX)%*%mX
  rho = 0.8
  
  ## Initial Global Mean and Global Covariance From MFVB
  vmu_adjust = params$vmu_til
  mSigma_adjust = params$mSigma_til
  
  # Later will be changed
  a_til = params$a_til
  b_til = params$b_til
  c = lambda * (exp(lgamma(a_til+0.5) - lgamma(a_til) - 0.5*log(b_til)))
  
  # Record the local parameter
  va = list()
  vb = list()
  vZ = c()
  
  # Initial value Local parameters
  vtheta = c(vmu_adjust)
  
  # Get unique pair of variable combination
  pairs = t(combn(unique(c(1:p)),2))
  
  
  for (ITER in 1:MAXITER)
  {
    ## Local Update
    for (j in 1:dim(pairs)[1])
    {
      # Store parameter for previous iteration
      vmu_old = vmu_adjust
      mSigma_old = mSigma_adjust
      curr_pair = c(pairs[j,1],pairs[j,2])
      
      # Define some constant
      mSigma_jj_inv = solve(mSigma_old[curr_pair,curr_pair])
      
      mt = mSigma_old[-curr_pair,curr_pair] %*% mSigma_jj_inv
      vs = as.matrix(vmu_old[-curr_pair]) - mt %*% vmu_old[curr_pair]
      
      
      # Local Update
      ## Update local parameter
      A = a_til/b_til*(XTX[curr_pair,curr_pair] + 0.5*XTX[curr_pair,-curr_pair]%*%mt + 0.5* t(mt)%*%XTX[-curr_pair,curr_pair])
      b = a_til/b_til*t(mX[,curr_pair])%*%(vy-mX[,-curr_pair]%*%vs)
      
      
      ## Calculate Local mean and variance
      vlocal_mean = emlasso(A,b,c)
      mlocal_var = vmlasso(A,b,c)
      
      # Record local parameter
      va[[j]] = A
      vb[[j]] = b
      vZ[j] = zmlasso(A,b,c)
      
      
      # Global Update
      ## Update Mean
      vmu_adjust[curr_pair] = vlocal_mean
      vmu_adjust[-curr_pair] = vmu_adjust[-curr_pair] + mSigma_adjust[-curr_pair,curr_pair] %*%  
        mSigma_jj_inv  %*% (vlocal_mean - vmu_old[curr_pair])
      
      
      ## Update Covariance
      mSigma_adjust[curr_pair,curr_pair] = mlocal_var
      mSigma_adjust[curr_pair,-curr_pair] = mlocal_var %*% mSigma_jj_inv %*%  mSigma_old[curr_pair,-curr_pair]
      mSigma_adjust[-curr_pair,curr_pair] = t(mSigma_adjust[curr_pair,-curr_pair])
      mSigma_adjust[-curr_pair,-curr_pair] =  as.matrix(mSigma_old[-curr_pair,-curr_pair]) + 
        mSigma_old[-curr_pair,curr_pair] %*% mSigma_jj_inv %*% 
        (mlocal_var - mSigma_old[curr_pair,curr_pair]) %*% mSigma_jj_inv %*%
        mSigma_old[curr_pair,-curr_pair]      
      
      
      
      # Damping
      # mSigma_adjust = rho* mSigma_adjust + (1-rho)*mSigma_old
      
      
    }
    
    # Update Theta
    vtheta_old <- vtheta
    vtheta <- c(vmu_adjust)
    
    ## Check stopping criterion    
    err <- max(abs(vtheta - vtheta_old))
    cat("ITER=",ITER,"err=",err,"\n")
    if (err < TOL) {
      break;
    }
  }
  # Return Local and Global parameter
  return(list("vmu_til" = vmu_adjust,"mSigma_til" = mSigma_adjust, "A" = va,"b" = vb,"c"= c ,"Z" = vZ))
  
}