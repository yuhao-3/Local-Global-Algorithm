lasso_gibbs_2B <- function(mX,vy,u,v, a = 1, b = 1.2, max_iter)
{
  n <- nrow(mX)
  p <- ncol(mX)
  
  # Storage for MCMC
  mBeta_mcmc   <- matrix(NA,max_iter,p)
  sigma2_mcmc  <- rep(NA, max_iter)
  lambda2_mcmc <- rep(NA, max_iter)
  
  mMode_mcmc   <- matrix(NA,max_iter,p)
  
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
    for(j in 1:p)
    {
      num <- XTy[j] - XTX[j,-j]%*%vbeta[-j]
      
      ## Note that XTX[j,j] = sum(x[,j]^2)
      a_val <- XTX[j,j]/sigma2
      c_val <- sqrt(lambda2/sigma2)
      b_val <- as.vector(num/sigma2)
      
      vbeta[j] <- rlasso_fast(a = a_val, b = b_val, c = c_val)
      
      mode <- mlasso(a_val, b = b_val, c = c_val)
      
      mMode_mcmc[iter,j] <- mode
      
      #beta <- arms(1, 
      #                function(x) -0.5*a_val*x*x + b_val*x - c_val*abs(x), 
      #                lower = -1000, upper = 1000, metropolis = TRUE)
      #vbeta[j] <- beta
      
    }
    vbeta <- matrix(vbeta)
    
    sigma <- sqrt(sigma2)
    u_til <- u + 0.5*p
    sum_abs_vbeta <- sum(abs(vbeta))
    lambda2 <- arms(1, 
                    function(x) (u_til - 1)*log(x) - v*x - sqrt(x)*sum_abs_vbeta/sigma, 
                    lower = 0, upper = 1000, metropolis = TRUE)
    lambda <- sqrt(lambda2)
    
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vbeta)^2) 
    sigma2 <- arms(1, 
                   function(x) -(a_til+1)*log(x) - b_til/x - lambda*sum_abs_vbeta/sqrt(x), 
                   lower = 0, upper = 1000, metropolis = TRUE)
    
    
    mBeta_mcmc[iter,] <- vbeta
    sigma2_mcmc[iter] <- sigma2
    lambda2_mcmc[iter] <- lambda2
    
    #if(iter%in%seq(0,max_iter,by=100)) print(iter)
  }
  
  
  output <- list()
  output$mBeta_mcmc = mBeta_mcmc
  output$sigma2_mcmc = sigma2_mcmc
  output$lambda2_mcmc = lambda2_mcmc
  output$mMode_mcmc = mMode_mcmc
  return(output)
}
