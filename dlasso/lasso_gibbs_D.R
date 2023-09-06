
source("sample_sigma2_and_lambda2.R")

################################################################################


lasso_gibbs_D <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter, sI, vbeta, lambda2, sigma2)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  sJ <- (1:p)[-sI]
  
  p1 <- length(sI)
  p0 <- length(sJ)
  
  mX1 <- mX[,sI]
  mX0 <- mX[,sJ]
  
  # Create storage for Rao-Blackwellisation
  mA_rw <- matrix(NA,max_iter,p0)
  mB_rw <- matrix(NA,max_iter,p0)
  mC_rw <- matrix(NA,max_iter,1)
  mM <- matrix(NA,max_iter,p1)
  mV <- matrix(NA,max_iter,p1)
  
  vy <- as.matrix(vy)
  XTy <- t(mX)%*%vy
  yTy <- t(vy)%*%vy
  
  #dgXTX <- matrix(1,1,n)%*%(mX*mX)
  #dgXTX1 <- dgXTX[sI]
  dgXTX0 <- matrix(1,1,n)%*%(mX0*mX0)
  X1Ty <- t(mX1)%*%vy
  X0Ty <- t(mX0)%*%vy
  X1TX1 <- t(mX1)%*%mX1
  X1TX0 <- t(mX1)%*%mX0
  #X0TX  <- t(mX0)%*%mX
  
  # Set the current values of the parameters 
  va    <- rep(1,p1)
 
  a_til <- a + 0.5*(n + p)
  u_til <- u + 0.5*p

  accept_sigma2 <- 0
  accept_lambda2 <- 0
  
  for (iter in 1:max_iter)
  {
    ############################################################################
    
    vbeta1 <- vbeta[sI]
    vbeta0 <- vbeta[sJ]
    
    vB <- lambda2*vbeta1*vbeta1/sigma2 
    vnu <- 1/sqrt(vB)
    for (j in 1:p1) {
      #va[j] <- rgig(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[i] <- rgig(1, lambda=0.5, chi=1, psi=vB[i])
      #vb[j] <- rinvgauss(1, mean=vmu[j])
      va[j] <- 1/statmod::rinvgauss(1, vnu[j], 1)
    }
    #vb <- as.vector(vb)
    va <- as.vector(va)
    
    ############################################################################
    
    # Sample from vbeta[sI]
    mQ_inv <- X1TX1 + diag(lambda2/va)
    mQ <- solve(mQ_inv,tol=1.0E-99)
    mSigma <- sigma2*mQ
    vmu1 <- mQ%*%X1Ty
    mP <- mQ%*%X1TX0
    vmu <- vmu1 - mP%*%vbeta0
    vmu <- as.vector(vmu)
    vbeta1 <- rmvnorm(n = 1, mean = vmu, sigma = mSigma)
    vbeta1 <- as.vector(vbeta1)
    
    mM[iter,] <- vmu
    mV[iter,] <- diag(mSigma)
    
    ############################################################################

    vu <- runif(p0)
    va_val <- (dgXTX0 - matrix(1,1,p1)%*%(X1TX0*mP))/sigma2
    vb_val <- rep(0,p0)
    c_val <- sqrt(lambda2/sigma2)
    
    mX0_hat <- mX0 - mX1%*%mP
 
    vy_hat <- mX1%*%vmu1 + mX0_hat%*%vbeta0
    
    for(j in 1:p0)
    {
      vy_hat_mj <- vy_hat - mX0_hat[,j]*vbeta0[j]
      vb_val[j] <- (X0Ty[j] - sum(mX0[,j]*vy_hat_mj))/sigma2
      vbeta0[j]  <- qlasso_fast(vu[j], va_val[j], vb_val[j], c_val)
      vy_hat <- vy_hat_mj + mX0_hat[,j]*vbeta0[j]
    }
    
    vy_hat <- mX1%*%vbeta1 + mX0%*%vbeta0
    
    #print(va_val)
    #print(vb_val)
    #print(vbeta0)

    # Store lasso parameters for Rao-Blackwellisation
    mA_rw[iter,] <- va_val
    mB_rw[iter,] <- vb_val
    mC_rw[iter]  <- c_val
    
    ############################################################################
    
    # Metropolis-Hastings step for lambda2
    
    sigma <- sqrt(sigma2)
    lambda <- sqrt(lambda2)
    sum_abs_vbeta <- sum(abs(vbeta0))
    RSS <- sum((vy - vy_hat)^2)
    sum_vbeta2_on_va <- sum(vbeta1*vbeta1/va)
    
    a_val <- u_til - 1
    b_val <- v + 0.5*sum_vbeta2_on_va/sigma2
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
    b_val <- b_til + 0.5*sum_vbeta2_on_va*lambda2
    c_val <- lambda*sum_abs_vbeta
    
    res <- sample_sigma2(sigma2, a_val, b_val, c_val) 
    sigma2 <- res$sigma2
    accept_sigma2 <- accept_sigma2 + res$accept
    
    ############################################################################
    
    vbeta[sI] <- vbeta1
    vbeta[sJ] <- vbeta0
    
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
