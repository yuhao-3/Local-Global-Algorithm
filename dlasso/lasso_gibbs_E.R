
source("sample_sigma2_and_lambda2.R")

################################################################################


lasso_gibbs_E <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter, sI, vbeta, lambda2, sigma2)
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
  mM <- matrix(NA,max_iter,p)
  mV <- matrix(NA,max_iter,p)
  va_til_rb   <- matrix(NA,max_iter,1)
  vb_til_rb   <- matrix(NA,max_iter,1)
  vu_til_rb   <- matrix(NA,max_iter,1)
  vv_til_rb   <- matrix(NA,max_iter,1)
  
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
  
  vmu_til <- rep(0,p0)
  vsigma2_til <- rep(1,p0)

  accept_sigma2 <- 0
  accept_lambda2 <- 0
  
  for (iter in 1:max_iter)
  {
    ############################################################################
    
    vbeta1 <- vbeta[sI]
    vbeta0 <- vbeta[sJ]
    
    vB <- lambda2*vbeta*vbeta/sigma2 
    vnu <- 1/sqrt(vB)
    for (j in 1:p) {
      #va[j] <- rgig(1, chi = vB[j], psi = 1, lambda = 0.5)
      #va[i] <- rgig(1, lambda=0.5, chi=1, psi=vB[i])
      #vb[j] <- rinvgauss(1, mean=vmu[j])
      va[j] <- 1/statmod::rinvgauss(1, vnu[j], 1)
    }
    #vb <- as.vector(vb)
    va <- as.vector(va)
    vd <- lambda2/va
    
    ############################################################################
    
    # Sample from vbeta[sI]
    mQ_inv <- X1TX1 + diag(vd[sI])
    mQ <- solve(mQ_inv,tol=1.0E-99)
    mSigma <- sigma2*mQ
    vmu1 <- mQ%*%X1Ty
    mP <- mQ%*%X1TX0
    vmu <- vmu1 - mP%*%vbeta0
    vmu <- as.vector(vmu)
    vbeta1 <- rmvnorm(n = 1, mean = vmu, sigma = mSigma)
    vbeta1 <- as.vector(vbeta1)
    
    mM[iter,sI] <- vmu
    mV[iter,sI] <- diag(mSigma)
    
    ############################################################################

    vu <- runif(p0)
    va_val <- 1/(dgXTX0 - matrix(1,1,p1)%*%(X1TX0*mP) + vd[sJ])
    vsigma_til <- sqrt(sigma2*va_val)
    
    mX0_hat <- mX0 - mX1%*%mP
    vy_hat <- mX1%*%vmu1 + mX0_hat%*%vbeta0
    
    for(j in 1:p0)
    {
      vy_hat_mj <- vy_hat - mX0_hat[,j]*vbeta0[j]
      vmu_til[j] <- (X0Ty[j] - sum(mX0[,j]*vy_hat_mj))*va_val[j]
      vbeta0[j] <- rnorm(1, mean = vmu_til[j], sd = vsigma_til[j])
      vy_hat <- vy_hat_mj + mX0_hat[,j]*vbeta0[j]
    }
    
    mM[iter,sJ] <- vmu_til
    mV[iter,sJ] <- vsigma2_til
    
    vbeta[sI] <- vbeta1
    vbeta[sJ] <- vbeta0

    ############################################################################
    
    # Metropolis-Hastings step for lambda2
    
    vy_hat <- mX%*%vbeta 
    sum_vbeta2_on_va <- sum(vbeta*vbeta/va)
  
    b_til <- b + 0.5*sum((vy - vy_hat)^2) + 0.5*lambda2*sum_vbeta2_on_va
    b_til <- as.vector(b_til)
    sigma2 <- rinvgamma(1, shape=a_til, rate=b_til)

    v_til <- v + 0.5*sum_vbeta2_on_va/sigma2
    lambda2 <- rgamma(1, shape=u_til, rate=v_til)
    
    ############################################################################
    
    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2
    
    va_til_rb[iter] <- a_til
    vb_til_rb[iter] <- b_til
    vu_til_rb[iter] <- u_til
    vv_til_rb[iter] <- v_til
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
      mM = mM,
      mV = mV
    )
  )
}


