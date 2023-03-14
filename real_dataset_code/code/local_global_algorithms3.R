###############################################################################
source(here("code","lasso_distribution","bivariate_lasso_distribution.R"))

###############################################################################

check_posdef <- function(mA) {
  val <- TRUE
  if (any(is.na(mA))) {
    val <- FALSE
    return(val)
  }
  if (any(!is.finite(mA))) {
    val <- FALSE
    return(val)
  }
  if (!all(eigen(mA)$values>0)) {
    val <- FALSE
    return(val)
  }  
  return(val)
}



# When bem coefficient is 0, then use univariate, otherwise, use bivariate for other pairs
local_global_algorithm_3 <- function(vy, mX, lambda, vmu_init, mSigma_init, a_til, b_til, damp, vbem_beta)
{
  ## Initialization
  MAXITER = 500
  TOL = 1.0E-8
  n = nrow(mX)
  p = ncol(mX)
  
  VERBOSE <- TRUE
  
  XTX = t(mX)%*%mX
  XTy = t(mX)%*%vy
  yTy = t(vy)%*%vy
  
  rho = 0.8
  
  ## Initial Global Mean and Global Covariance From MFVB
  vmu_glob <- vmu_init
  mSigma_glob <- mSigma_init
  
  # The "c" parameter in the bivariate lasso is constant for the time being.
  c_val = lambda*(exp(lgamma(a_til+0.5) - lgamma(a_til) - 0.5*log(b_til)))
  #cat("c_val=",c_val,"\n")
  
  E_inv_sigma2 <- a_til/b_til
  #cat("E_inv_sigma2=",E_inv_sigma2,"\n")
  
  # Record the local parameter
  lmA = list()
  lvb = list()
  
  
  # Initial value Local parameters
  vtheta = c(vmu_glob)
  
  # Get index of variable with zero coefficient
  zero_var = which(vbem_beta == 0)
  
  # Get unique pair of non-zero variable combination
  mPairs = t(combn(unique(setdiff(c(1:p),zero_var)),2))
  
  
  for (ITER in 1:MAXITER)
  {
    # ## Local Update
    # for (j in 1:length(zero_var))
    # {
    #   # Store parameter for previous iteration
    #   vmu_old = vmu_glob
    #   mSigma_old = mSigma_glob
    #   
    #   # Define some constant
    #   j_len = length(vmu_glob[j])
    #   mSigma_jj_inv = matrix(solve(mSigma_glob[j,j]),j_len,j_len)
    #   mt = matrix(mSigma_glob[-j,j]) %*% mSigma_jj_inv
    #   vs = matrix(vmu_glob[-j]) - mt %*% matrix(vmu_glob[j])
    #   
    #   # Local Update
    #   ## Update local parameter
    #   a = a_til/b_til*(matrix(XTX[j,j]) + t(XTX[j,-j]%*%mt))
    #   b = a_til/b_til*t(mX[,j])%*%(vy-mX[,-j]%*%vs)
    #   
    #   
    #   ## Calculate Local mean and variance
    #   vlocal_mean = elasso(a,b,c)
    #   mlocal_var = vlasso(a,b,c)
    #   
    #   
    #   # Record local parameter
    #   lma[j] = a
    #   lvb[j] = b
    #   # vZ[j] = zlasso(a,b,c)
    #   
    #   # Global Update
    #   ## Update Mean
    #   
    #   vmu_glob[j] = vlocal_mean
    #   vmu_glob[-j] = matrix(vmu_glob[-j]) + matrix(mSigma_glob[-j,j]) %*%  mSigma_jj_inv  %*% matrix(vlocal_mean - vmu_old[j])
    #   
    #   ## Update Covariance
    #   mSigma_glob[j,j] = mlocal_var
    #   mSigma_glob[j,-j] = mlocal_var %*% mSigma_jj_inv %*%  t(matrix(mSigma_old[j,-j]))
    #   mSigma_glob[-j,j] = t(mSigma_glob[j,-j])
    #   mSigma_glob[-j,-j] = matrix(mSigma_old[-j,-j],p-j_len,p-j_len) + matrix(mSigma_old[-j,j]) %*% mSigma_jj_inv %*% (mlocal_var - matrix(mSigma_old[j,j]))  %*%mSigma_jj_inv %*% t(matrix(mSigma_old[j,-j]))
    #   
    # }
    
    
    
    for (j in 1:nrow(mPairs))
    {
      vmu_old <- vmu_glob
      mSigma_old <- mSigma_glob
      
      pair <- mPairs[j,]
      
      
      
      ##########################################################################
      
      # Local update
      
      PASS <- TRUE
      
      if (!check_posdef(mSigma_glob[pair,pair])) {
        print("Something has gone terribly wrong 1")
        print(mSigma_glob[pair,pair])
        stop()
        PASS <- FALSE
      }
      
      mOmega <- solve(mSigma_glob[pair,pair], tol=1.0E-99)
      
      if (!check_posdef(mOmega)) {
        print("Something has gone terribly wrong 2")
        print(mOmega)
        stop()
        PASS <- FALSE
      }
      
      
      mT = mSigma_glob[-pair,pair]%*%mOmega
      vs = vmu_glob[-pair] - mT%*%vmu_glob[pair]
      
      # Calculate the parameters of the lasso distribution
      mA_val = E_inv_sigma2*(XTX[pair,pair] + XTX[pair,-pair]%*%mT)
      vb_val = E_inv_sigma2*(XTy[pair] - XTX[pair,-pair]%*%vs)
      
      # Make mA_val symmetric
      mA_val <- (mA_val + t(mA_val))/2 
      
      if (!check_posdef(mA_val)) {
        print("Something has gone terribly wrong 3")
        print(mA_val)
        stop()
        PASS <- FALSE
      }
      
      
      
      
      #
      
      if (PASS) {
        
        # Calculate normalizing constant of the bivariate lasso distribution
        #log_Z <- zbvlasso(mA_val,vb_val,c_val, logarithm=TRUE) 
        
        # Calculate the mean and covariance of the bivariate lasso distribution
        vmu_star    <- ebvlasso(mA_val,vb_val,c_val) 
        mSigma_star <- vbvlasso(mA_val,vb_val,c_val) 
        
        #if (any(abs(vmu_star)>0.6)) 
        #{
        #  PASS <- FALSE
        #}
        
        if (!check_posdef(mSigma_star)) {
          PASS <- FALSE
        }
        
        if (VERBOSE) {
          
          
          cat("pair=",pair,"\n")
          cat("mA_val=",mA_val,"\n")
          cat("vb_val=",vb_val,"\n")
          cat("c_val=",c_val,"\n")
          cat("mSigma_star=",mSigma_star,"\n")
          
          print(sbvlasso(mA_val,vb_val,c_val))
          
          print(ebvlasso(mA_val,vb_val,c_val,verbose=TRUE))
          print(vbvlasso(mA_val,vb_val,c_val,verbose=TRUE))
          
          
        }
        
        if (!PASS) {
          print("Something has gone terribly wrong 4")
          stop()
        }
        
        #cat("mSigma_star=",mSigma_star,"\n")
        if (PASS) 
        {
          # Apply dampening so things don't change quickly
          vmu_star    <- damp*vmu_star + (1-damp)*vmu_old[pair] 
          mSigma_star <- damp*mSigma_star + (1-damp)*mSigma_old[pair,pair]  
          
          
          # Record local parameter
          lmA[[j]] = mA_val
          lvb[[j]] = vb_val
          #vZ[j] = log_Z
          
          ##########################################################################
          
          # Global update
          
          ## Update Mean
          vmu_glob[pair] = vmu_star
          vmu_glob[-pair] = vmu_glob[-pair] + mSigma_glob[-pair,pair]%*%mOmega%*%(vmu_star - vmu_glob[pair])
          
          ## Update Covariance
          mSigma_glob[pair,pair] = mSigma_star
          mSigma_glob[pair,-pair] = mSigma_star%*%mOmega%*% mSigma_old[pair,-pair]
          mSigma_glob[-pair,pair] = t(mSigma_glob[pair,-pair])
          
          mSigma_glob[-pair,-pair] = mSigma_glob[-pair,-pair] + 
            mSigma_old[-pair,pair]%*%mOmega%*%(mSigma_star - mSigma_old[pair,pair])%*%mOmega%*%mSigma_old[pair,-pair]      
          
          if (!check_posdef(mSigma_glob)) {
            print("Something has gone terribly wrong 5")
            print(mSigma_glob)
            stop()
          }
        } 
      }
      ##########################################################################
    }
    
    # Update Theta
    vtheta_old <- vtheta
    vtheta <- c(vmu_glob)
    
    ## Check stopping criterion    
    err <- max(abs(vtheta - vtheta_old))
    cat("ITER=",ITER,"err=",err,"\n")
    if (err < TOL) {
      break;
    }
    
    
  }
  # Return Local and Global parameter
  return(list("vmu_til" = vmu_glob, "mSigma_til" = mSigma_glob, "lmA" = lmA,"lvb" = lvb,"c_val"= c_val, mPairs=mPairs))
  
}







