lasso_gibbs_2A <- function(mX,vy,a, b,u,v, max_iter, vbeta, lambda2, sigma2)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  va_til_rb   <- matrix(NA,max_iter,1)
  vb_til_rb   <- matrix(NA,max_iter,1)
  vu_til_rb   <- matrix(NA,max_iter,1)
  vv_til_rb   <- matrix(NA,max_iter,1)
  mM_rb <- matrix(NA,max_iter,p)
  mV_rb <- matrix(NA,max_iter,p)
  
  # Not needed since we are not doing inference on va
  #mA_mcmc     <- array(0, dim = c(p,iter))
  
  vy <- as.matrix(vy)
  
  if (n>p) {
    XTX <- t(mX)%*%mX
  } else {
    dgXTX <- matrix(1,1,n)%*%(mX*mX)
  }
  
  YTX <- t(vy)%*%mX
  XTy <- t(mX)%*%vy
  
  # Set the current values of the parameters 
  va    <- rep(1,p)
   
  
  for (iter in 1:max_iter)
  {
    sigma <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    
    vnu <- sigma/(lambda*abs(vbeta))  
    #vB <- lambda2*vbeta*vbeta/sigma2 
    va <-  1/rinvgauss(p, mean=vnu, rep(1,p))
    
    
    #for (j in 1:p) {
    #  va[j] <- statmod::rinvgauss(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[j] <- rgig(1, lambda=0.5, chi=1, psi=vB[j])
      #vb[j] <- rinvgauss(1, mean=vmu[j])
      #vb[j] <- rinvGauss(1, nu=vnu[j], lambda=1)
    #}
    #vb <- as.vector(vb)
    va <- as.vector(va)
    
    if (n>p) {
      for(j in 1:p){
        num <- XTy[j] - XTX[j,-j]%*%vbeta[-j]
        denom <- XTX[j,j] + lambda2/va[j]
        mu_til <- num/denom
        sigma2_til <- sigma2/denom
        vbeta[j] <- rnorm(1, mean = mu_til, sd = sqrt(sigma2_til))
        
        mM_rb[iter,j] <- mu_til
        mV_rb[iter,j] <- sigma2_til
      }
    } else {
      vy_hat <- mX%*%vbeta
      for(j in 1:p){
        vy_hat_mj <- vy_hat - mX[,j]*vbeta[j]
        num <- XTy[j] - sum(mX[,j]*vy_hat_mj)
        denom <- dgXTX[j] + lambda2/va[j]
        mu_til <- num/denom
        sigma2_til <- sigma2/denom
        vbeta[j] <- rnorm(1, mean = mu_til, sd = sqrt(sigma2_til))
        vy_hat <- vy_hat_mj + mX[,j]*vbeta[j]
        
        mM_rb[iter,j] <- mu_til
        mV_rb[iter,j] <- sigma2_til
      }
    }
    vbeta <- matrix(vbeta)
    
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vbeta)^2) + 0.5*lambda2*sum(vbeta*vbeta/va)
    b_til <- as.vector(b_til)
    
    sigma2 <- rinvgamma(1, shape=a_til, rate=b_til)
    #sigma <- sqrt(sigma2)
    

    
    u_til <- u + 0.5*p
    v_til <- v + 0.5*sum(vbeta*vbeta/va)/sigma2
    lambda2 <- rgamma(1, shape=u_til, rate=v_til)

    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2
    
    va_til_rb[iter] <- a_til
    vb_til_rb[iter] <- b_til
    vu_til_rb[iter] <- u_til
    vv_til_rb[iter] <- v_til

    
    #if(iter%in%seq(0,max_iter,by=100)) print(iter)
  }
  
  return(
    list(
      mBeta = mBeta_mcmc,
      vsigma2 = sigma2_mcmc,
      vlambda2 = lambda2_mcmc,
      va_til = va_til_rb,
      vb_til = vb_til_rb,
      vu_til = vu_til_rb,
      vv_til = vv_til_rb,
      mM = mM_rb,
      mV = mV_rb
    )
  )
}
  
