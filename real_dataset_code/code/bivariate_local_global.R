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


###############################################################################

# Bivariate Local and Global Iterations 
# (Update each pair of variables at a time) 

bivariate_local_global <- function(vy, mX, lambda, vmu_init, mSigma_init, a_til, b_til)
{
  ## Initialization
  MAXITER = 500
  TOL = 1.0E-8
  n = nrow(mX)
  p = ncol(mX)
  
  XTX = t(mX)%*%mX
  XTy = t(mX)%*%vy
  yTy = t(vy)%*%vy
  
  rho = 0.8
  
  ## Initial Global Mean and Global Covariance From MFVB
  vmu_glob <- vmu_init
  mSigma_glob <- mSigma_init
  
  # The "c" parameter in the bivariate lasso is constant for the time being.
  c_val = lambda*(exp(lgamma(a_til+0.5) - lgamma(a_til) - 0.5*log(b_til)))
  cat("c_val=",c_val,"\n")
  
  E_inv_sigma2 <- a_til/b_til
  cat("E_inv_sigma2=",E_inv_sigma2,"\n")
  
  # Record the local parameter
  va = list()
  vb = list()
  vZ = c()
  
  # Initial value Local parameters
  vtheta = c(vmu_glob)
  
  # Get unique pair of variable combination
  mPairs = t(combn(unique(c(1:p)),2))
  
  
  for (ITER in 1:MAXITER)
  {
    ## Local Update
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
        PASS <- FALSE
      }
      
      mOmega <- solve(mSigma_glob[pair,pair], tol=1.0E-99)
      
      if (!check_posdef(mOmega)) {
        print("Something has gone terribly wrong 2")
        print(mOmega)
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
        PASS <- FALSE
      }
      
    

      
      #
      
      if (PASS) {
      
        # Calculate normalizing constant of the bivariate lasso distribution
        #log_Z <- zbvlasso(mA_val,vb_val,c_val, logarithm=TRUE) 
        
        # Calculate the mean and covariance of the bivariate lasso distribution
        vmu_star    <- ebvlasso(mA_val,vb_val,c_val) 
        mSigma_star <- vbvlasso(mA_val,vb_val,c_val) 
        
        if (!check_posdef(mSigma_star)) {
          print("Something has gone terribly wrong 4")
          
          cat("pair=",pair,"\n")
          cat("mA_val=",mA_val,"\n")
          cat("vb_val=",vb_val,"\n")
          cat("c_val=",c_val,"\n")
          cat("mSigma_star=",mSigma_star,"\n")
          
          print(sbvlasso(mA_val,vb_val,c_val))
          
          
          PASS <- FALSE
        }
        
        #cat("mSigma_star=",mSigma_star,"\n")
        if (PASS) 
        {
          
          # Apply dampening so things don't change quickly
          alpha <- 0.75
          vmu_star    <- alpha*vmu_star + (1-alpha)*vmu_old[pair] 
          mSigma_star <- alpha*mSigma_star + (1-alpha)*mSigma_old[pair,pair]  
          
          
          # Record local parameter
          va[[j]] = mA_val
          vb[[j]] = vb_val
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
            break;
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
  return(list("vmu_til" = vmu_glob, "mSigma_til" = mSigma_glob, "A" = va,"b" = vb,"c"= c ,"Z" = vZ))
  
}