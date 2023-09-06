source("sample_sigma2_and_lambda2.R")

################################################################################


lasso_gibbs_2B <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter, vbeta_init, lambda2_init, sigma2_init)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  # Store statistics for Rao-Blackwellisation
  mA_rw <- matrix(NA,max_iter,p)
  mB_rw <- matrix(NA,max_iter,p)
  mC_rw <- matrix(NA,max_iter,p)
  
  # Calculate summary statistics
  vy <- as.matrix(vy)
  if (n>p) {
    XTX <- t(mX)%*%mX
  } 
  yTy <- t(vy)%*%vy
  XTy <- t(mX)%*%vy
  dgXTX <- matrix(1,1,n)%*%(mX*mX)
  
  a_til <- a + 0.5*(n + p)
  u_til <- u + 0.5*p
  
  accept_sigma2 <- 0
  accept_lambda2 <- 0
  
  vbeta = vbeta_init
  sigma2 = sigma2_init
  lambda2 = lambda2_init
  
  for (iter in 1:max_iter)
  {
    va_val <- dgXTX/sigma2
    c_val  <- sqrt(lambda2/sigma2)
    
    mA_rw[iter,] <- va_val
    mC_rw[iter,] <- c_val
    
    vu <- runif(p)
    
    if (n>p) {
      for(j in 1:p)
      {
        num <- XTy[j] - XTX[j,-j]%*%vbeta[-j]
        b_val <- num/sigma2
        vbeta[j] <- qlasso_fast(vu[j],va_val[j], b = b_val, c = c_val) #rlasso_fast(a = va_val[j], b = b_val, c = c_val)
        mB_rw[iter,j] <- b_val
      }
    } else {
      vy_hat <- mX%*%vbeta
      for(j in 1:p)
      {
        vy_hat_mj <- vy_hat - mX[,j]*vbeta[j]
        num <- XTy[j] - sum(mX[,j]*vy_hat_mj)
        b_val <- num/sigma2
        vbeta[j] <- qlasso_fast(vu[j],va_val[j], b = b_val, c = c_val)  #rlasso_fast(a = va_val[j], b = b_val, c = c_val)
        vy_hat <- vy_hat_mj + mX[,j]*vbeta[j]
        mB_rw[iter,j] <- b_val
      }
    }
    vbeta <- matrix(vbeta)
    
    ############################################################################
    
    sigma <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    sum_abs_vbeta <- sum(abs(vbeta))
    RSS <- sum((vy - mX%*%vbeta)^2)
    
    a_val <- u_til - 1
    b_val <- v
    c_val <- sum_abs_vbeta/sigma
    
    # Slice sampler for lambda2
    lambda2 <- slice_sample(lambda2 , a_val, b_val, c_val)
    
    # Metropolis-Hastings step for lambda2
    #res <- sample_lambda2(lambda2, a_val, b_val, c_val) 
    #lambda2 <- res$lambda2
    #accept_lambda2 <- accept_lambda2 + res$accept
    
    ############################################################################
    
    
    
    sigma  <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    
    a_val <- (a_til-1)
    b_val <- b + 0.5*RSS
    c_val <- lambda*sum_abs_vbeta
    
    # Slice sampler for tau and then invert.
    tau <- 1/sigma2
    tau <- slice_sampler(tau, a_val, b_val, c_val)
    sigma2 <- 1/tau
    
    # Metropolis-Hastings step for sigma2
    #res <- sample_sigma2(sigma2, a_val, b_val, c_val) 
    #sigma2 <- res$sigma2
    #accept_sigma2 <- accept_sigma2 + res$accept
 
    ############################################################################
    
    if ((iter%%1000)==0)
      cat("iter=",iter, "sigma2=",sigma2, "lambda2=", lambda2, "\n")
    
    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2
  }
  
  output <- list()
  output$mBeta_mcmc = mBeta_mcmc
  output$sigma2_mcmc = sigma2_mcmc
  output$lambda2_mcmc = lambda2_mcmc
  output$accept_sigma2 = accept_sigma2
  output$accept_lambda2 = accept_lambda2
  output$mA_rw = mA_rw
  output$mB_rw = mB_rw
  output$mC_rw = mC_rw
  
 
  return(output)
}




density_rw_beta2 <- function(x, mA, mB, mC,j) 
{
  y <- x*0
  N <- nrow(mA)
  for (i in 1:N) {
    y <- y + dlasso(x, a = mA[i,j], b = mB[i,j], c = mC[i,j])
  }
  y <- y/N
  return(list(x=x,y=y))
}
