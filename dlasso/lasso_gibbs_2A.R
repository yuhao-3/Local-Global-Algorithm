lasso_gibbs_2A <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  # Not needed since we are not doing inference on va
  #mA_mcmc     <- array(0, dim = c(p,iter))
  
  vy <- as.matrix(vy)
  
  XTX <- t(mX)%*%mX
  YTX <- t(vy)%*%mX
  XTy <- t(mX)%*%vy
  
  # Set the current values of the parameters 
  va    <- rep(1,p)
  vbeta <- matrix(0,p,1)
  lambda2 <- 14
  sigma2 <- 1
  
  for (iter in 1:max_iter)
  {
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vbeta)^2) + 0.5*lambda2*sum(vbeta*vbeta/va)
    b_til <- as.vector(b_til)
    sigma2 <- rinvgamma(1, shape=a_til, rate=b_til)
    
    for(j in 1:p){
      num <- XTy[j] - XTX[j,-j]%*%vbeta[-j]
      denom <- XTX[j,j] + lambda2/va[j]
      mu_til <- num/denom
      sigma2_til <- sigma2/denom
      vbeta[j] <- rnorm(1, mean = mu_til, sd = sqrt(sigma2_til))
    }
    vbeta <- matrix(vbeta)
    
    for (j in 1:p) {
      va[j] <- rgig(1,1,lambda2*(vbeta[j])^2/sigma2,1/2)
    }
    
    u_til <- u + 0.5*p
    v_til <- v + 0.5*sum(vbeta*vbeta/va)/sigma2
    lambda2 <- rgamma(1,shape = u_til,rate = v_til)
    
    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2
    
    #if(iter%in%seq(0,max_iter,by=100)) print(iter)
  }
  
  
  output <- list()
  output$mBeta_mcmc = mBeta_mcmc
  output$sigma2_mcmc = sigma2_mcmc
  output$lambda2_mcmc = lambda2_mcmc
  return(output)
}
  
