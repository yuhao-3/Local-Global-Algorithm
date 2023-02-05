
f0_fun <- function(mu,sigma) {
  mu_til <- mu/sigma
  val <- mu*(pnorm(mu_til) - pnorm(-mu_til)) + 2*sigma*dnorm(mu_til)
  return(val)
}

f1_fun <- function(mu,sigma) {
  mu_til <- mu/sigma
  val <- (pnorm(mu_til) - pnorm(-mu_til))  
  return(val)
}

f2_fun <- function(mu,sigma) {
  mu_til <- mu/sigma
  val <-  2*dnorm(mu_til)/sigma
  return(val)
}

trapint <- function(xgrid, fgrid) 
{
  ng <- length(xgrid)
  xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
  fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
  integ <- sum(xvec * fvec)/2
  return(integ)
}

bayesian_lasso_spvb <- function(vy, mX, lambda, vmu_til, mSigma_til, sigma2_hat) 
{
  MAXITER <- 500
  TOL <- 1.0E-8
  
  a <- 0.0001
  b <- 0.0001
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  XTX <- t(mX)%*%mX 
  XTy <- t(mX)%*%vy
  
 
  vsigma2_til <- diag(mSigma_til)
  vsigma_til <- sqrt(vsigma2_til)
  
  E_sigma2inv <- 1/sigma2_hat
  E_sigmainv <- sqrt(E_sigma2inv)
  
 
  vtheta <- c(vmu_til)
  for (ITER in 1:MAXITER) 
  {
    # M-step
    c1 <- E_sigma2inv
    c2 <- lambda*E_sigmainv
    
    vd <- f2_fun(vmu_til,vsigma_til)
    vd <- as.vector(vd)
    mD <- diag(vd)
    
    mSigma_til  <- solve(c1*XTX + c2*mD, tol=1.0E-99)
    vsigma2_til <- diag(mSigma_til)
    vsigma_til  <- sqrt(vsigma2_til)
    
    vg <- c1*(XTy - XTX%*%vmu_til) - c2*f1_fun(vmu_til,vsigma_til) 
    vmu_til <- vmu_til + mSigma_til%*%vg
    

    a_til <- a + 0.5*(n + p)
    b_til <- lambda*sum(f0_fun(vmu_til,vsigma_til))
    c_til <- b + 0.5*sum((vy - mX%*%vmu_til)^2) + 0.5*sum(diag(XTX%*%mSigma_til))
     
    
    sigma2g <- seq(0.0001,1,length=100000)
    fg <- -(a_til+1)*log(sigma2g) - b_til/sqrt(sigma2g) - c_til/sigma2g 
    fg <- fg - max(fg)
    
    #plot(sigma2g,exp(fg),type="l")
    #ans <- readline()
    
     
    
    Z <- trapint(sigma2g,exp(fg))
    
    E_sigma2inv <- trapint(sigma2g,(1/sigma2g)*exp(fg))/Z
    E_sigmainv  <- trapint(sigma2g,(1/sqrt(sigma2g))*exp(fg))/Z
    
    vtheta_old <- vtheta
    vtheta <- vtheta_old + 1*c(vmu_til-vtheta_old)
    
    
    
    err <- max(abs(vtheta - vtheta_old))
    cat("ITER=",ITER,"err=",err,"E_sigma2inv=",E_sigma2inv,"\n")  
    if (err < TOL) {
      break;
    }
  }
  
  
  
  return(list(vmu_til=vmu_til, mSigma_til=mSigma_til, a_til=a_til, b_til=b_til, vd=vd, E_sigma2inv=E_sigma2inv, E_sigmainv=E_sigmainv))
}