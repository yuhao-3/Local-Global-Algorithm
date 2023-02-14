source(here("code","lasso_distribution","myLasso.R")) 
source(here("code","lasso_distribution","Multi-Lasso.R")) 
library(gtools)

## Local and Global Algorithm One(Update one variable at a time)
## If there is a scalar output, change it to 1*1 matrix output, otherwise remain same
## Expect an quick but less performance

local_global_algorithm_1 <- function(vy, mX, lambda, params)
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
  c = lambda * gamma(a_til+0.5)/(gamma(a_til)*sqrt(b_til))

  
  # Record the local parameter
  va = c()
  vb = c()
  vZ = c()
  
  # Initial value Local parameters
  vtheta <- c(vmu_adjust)
  
  for (ITER in 1:MAXITER)
  {
    ## Local Update
    for (j in 1:p)
    {
      # Store parameter for previous iteration
      vmu_old = vmu_adjust
      mSigma_old = mSigma_adjust
      
      # Define some constant
      j_len = length(vmu_adjust[j])
      mSigma_jj_inv = matrix(solve(mSigma_adjust[j,j]),j_len,j_len)
      mt = matrix(mSigma_adjust[-j,j]) %*% mSigma_jj_inv
      vs = matrix(vmu_adjust[-j]) - mt %*% matrix(vmu_adjust[j])
      
      # Local Update
      ## Update local parameter
      a = a_til/b_til*(matrix(XTX[j,j]) + t(XTX[j,-j]%*%mt))
      b = a_til/b_til*t(mX[,j])%*%(vy-mX[,-j]%*%vs)
      
      
      ## Calculate Local mean and variance
      vlocal_mean = elasso(a,b,c)
      mlocal_var = vlasso(a,b,c)
      
      
      # Record local parameter
      va[j] = a
      vb[j] = b
      vZ[j] = zlasso(a,b,c)
      
      # Global Update
      ## Update Mean
      
      vmu_adjust[j] = vlocal_mean
      vmu_adjust[-j] = matrix(vmu_adjust[-j]) + matrix(mSigma_adjust[-j,j]) %*%  mSigma_jj_inv  %*% matrix(vlocal_mean - vmu_old[j])
      
      ## Update Covariance
      mSigma_adjust[j,j] = mlocal_var
      mSigma_adjust[j,-j] = mlocal_var %*% mSigma_jj_inv %*%  t(matrix(mSigma_old[j,-j]))
      mSigma_adjust[-j,j] = t(mSigma_adjust[j,-j])
      mSigma_adjust[-j,-j] = matrix(mSigma_old[-j,-j],p-j_len,p-j_len) + matrix(mSigma_old[-j,j]) %*% mSigma_jj_inv %*% (mlocal_var - matrix(mSigma_old[j,j]))  %*%mSigma_jj_inv %*% t(matrix(mSigma_old[j,-j]))
      
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
  return(list("vmu_til" = vmu_adjust,"mSigma_til" = mSigma_adjust, "a" = va,"b" = vb,"c"= c ,"Z" = vZ))
  
}



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
  c = lambda * gamma(a_til+0.5)/(gamma(a_til)*sqrt(b_til))
  
  # Record the local parameter
  va = list()
  vb = list()
  vZ = c()
  
  # Initial value Local parameters
  vtheta <- c(vmu_adjust)
  
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
      A = a_til/b_til*(XTX[curr_pair,curr_pair] + XTX[curr_pair,-curr_pair]%*%mt)
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
      
      
      
      
      # Positive definiteness of local variance
      
      print(det(mlocal_var))
      print(det(mlocal_var - mSigma_old[curr_pair,curr_pair]))
      print(det(mSigma_adjust[-curr_pair,-curr_pair]))
      
      
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


