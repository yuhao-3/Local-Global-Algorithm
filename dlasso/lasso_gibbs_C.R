
source("sample_sigma2_and_lambda2.R")

################################################################################


lasso_gibbs_C <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter, sI, vbeta, lambda2, sigma2)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  mMode_mcmc   <- matrix(NA,max_iter,p)
  
  mA_rw <- matrix(NA,max_iter,p)
  mB_rw <- matrix(NA,max_iter,p)
  mC_rw <- matrix(NA,max_iter,p)
  
  mM <- matrix(NA,max_iter,length(sI))
  mV <- matrix(NA,max_iter,length(sI))
  
  vy <- as.matrix(vy)
  if (n>p) {
    XTX <- t(mX)%*%mX
  } 
  
  XTy <- t(mX)%*%vy
  yTy <- t(vy)%*%vy
  
  
  mX1 <- mX[,sI]
  mX0 <- mX[,-sI]
  
  dgXTX <- matrix(1,1,n)%*%(mX*mX)
  X1TX1 <- t(mX1)%*%mX1
  X1TX0 <- t(mX1)%*%mX0
  
  # Set the current values of the parameters 
  va    <- rep(1,length(sI))
 
  a_til <- a + 0.5*(n + p)
  u_til <- u + 0.5*p
  
  accept_sigma2 <- 0
  accept_lambda2 <- 0
  
  for (iter in 1:max_iter)
  {
    ############################################################################
    
    vB <- lambda2*vbeta[sI]*vbeta[sI]/sigma2
    vnu <- 1/sqrt(vB)
    for (j in 1:length(sI)) {
      #va[j] <- rgig(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[i] <- rgig(1, lambda=0.5, chi=1, psi=vB[i])
      va[j] <- 1/statmod::rinvgauss(1, vnu[j], 1)
      #vb[j] <- rinvGauss(1, nu=vnu[j], lambda=1)
    }
    #vb <- as.vector(vb)
    va <- as.vector(va)
    
    ############################################################################
    
    # Sample from vbeta[sI]
    mQ_inv <- X1TX1 + diag(lambda2/va)
    mQ <- solve(mQ_inv,tol=1.0E-99)
    mSigma <- sigma2*mQ
    vmu <- mQ%*%(XTy[sI] - X1TX0%*%vbeta[-sI])
    vmu <- as.vector(vmu)
    vbeta[sI] <- rmvnorm(n = 1, mean = vmu, sigma = mSigma)
    vbeta <- as.vector(vbeta)
    
    mM[iter,] <- vmu
    mV[iter,] <- diag(mSigma)
    
    ############################################################################

    vu <- runif(p)
    va_val <- dgXTX/sigma2
    vb_val <- rep(0,p)
    c_val <- sqrt(lambda2/sigma2)
    
    sJ <- (1:p)[-sI]
    
    if (n>p) {
      for(J in 1:length(sJ))
      {
        j <- sJ[J]
        num <- XTy[j] - XTX[j,-j]%*%vbeta[-j]
        vb_val[j] <- num/sigma2
        vbeta[j] <- qlasso_fast(vu[j], a = va_val[j], b = vb_val[j], c = c_val)
      }
    } else {
      vy_hat <- mX%*%vbeta
      for(J in 1:length(sJ))
      {
        j <- sJ[J]
        vy_hat_mj <- vy_hat - mX[,j]*vbeta[j]
        num <- XTy[j] - sum(mX[,j]*vy_hat_mj)
        vb_val[j] <- num/sigma2
        vbeta[j]  <- qlasso_fast(vu[j], a = va_val[j], b = vb_val[j], c = c_val)
        vy_hat <- vy_hat_mj + mX[,j]*vbeta[j]
      }
    }
    vbeta <- matrix(vbeta)
    
    # Store lasso parameters for Rao-Blackwellisation
    mA_rw[iter,] <- va_val
    mB_rw[iter,] <- vb_val
    mC_rw[iter,] <- c_val
    
    ############################################################################
    
    # Metropolis-Hastings step for lambda2
    
    sigma <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    sum_abs_vbeta <- sum(abs(vbeta[-sI]))
    RSS <- sum((vy - mX%*%vbeta)^2)
    sum_vbeta_sq_on_va <- sum(vbeta[sI]*vbeta[sI]/va)
    

    
    a_val <- u_til - 1
    b_val <- v + 0.5*sum_vbeta_sq_on_va/sigma2
    c_val <- sum_abs_vbeta/sigma
    

    
    res <- sample_lambda2(lambda2, a_val, b_val, c_val) 
    lambda2 <- res$lambda2
    accept_lambda2 <- accept_lambda2 + res$accept
    

    ############################################################################
    
    # Metropolis-Hastings step for sigma2
    
    sigma  <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    
    b_til <- b + 0.5*RSS
    
    a_val <- (a_til+1)
    b_val <- b_til + 0.5*sum_vbeta_sq_on_va*lambda2
    c_val <- lambda*sum_abs_vbeta
    
    res <- sample_sigma2(sigma2, a_val, b_val, c_val) 
    sigma2 <- res$sigma2
    accept_sigma2 <- accept_sigma2 + res$accept
 
    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2

    #if((iter%%100)==0) {
    #  print(iter)
    #  par(mfrow=c(1,2))
    #  plot(density(sigma2_mcmc[1:iter]))
    #  plot(density(lambda2_mcmc[1:iter]))
    #}
  }
  
  
  output <- list()
  output$mBeta_mcmc = mBeta_mcmc
  output$sigma2_mcmc = sigma2_mcmc
  output$lambda2_mcmc = lambda2_mcmc
  output$mMode_mcmc = mMode_mcmc
  output$accept_sigma2 = accept_sigma2
  output$accept_lambda2 = accept_lambda2
  output$mA_rw = mA_rw
  output$mB_rw = mB_rw
  output$mC_rw = mC_rw
  output$mM = mM
  output$mV = mV
  
 
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
