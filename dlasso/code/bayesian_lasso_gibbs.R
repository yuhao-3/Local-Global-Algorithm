library(mvtnorm)
library(invgamma)

# library(statmod)
# library(actuar)
# library(LaplacesDemon)
# library(rmutil)
# library(SuppDists)
# library(GeneralizedHyperbolic)
# library(GIGrvg)
# library(ghyp)

source("rinvgaussian.R")

################################################################################

bayesian_lasso_gibbs <-
  function(vy,
           mX,
           lambda,
           sigma2,
           vbeta_init,
           va,
           nsamples = 10000,
           a,
           b,
           trunc = 0
           # aux_sampler = c(
           #   "statmod",
           #   "actuar",
           #   "LaplacesDemon",
           #   "rmutil_rginvgauss",
           #   "rmutil_rinvgauss",
           #   "SuppDists",
           #   "GeneralizedHyperbolic",
           #   "GIGrvg",
           #   "ghyp"
           # )
           )
  {
    
  MAXITER <- nsamples
  
  if (is.numeric(lambda)) {
    lambda2 <- lambda*lambda
    USEPRIOR <- FALSE
  } else {
    if (is.list(lambda)) {
      USEPRIOR <- TRUE
      lambda_list <- lambda 
      lambda <- lambda_list$lambda_start
      lambda2 <- lambda*lambda
    }
  }
  
  n <- nrow(mX)
  p <- ncol(mX)
  
  if (p<n) {
    XTX <- t(mX)%*%mX 
  } else {
    mI_n <- diag(rep(1, n))
  }
  XTy <- t(mX)%*%vy
  yTy <- t(vy)%*%vy
  
  #vb <- rep(1,p)
  va <- rep(1,p)
  vbeta <- vbeta_init
  
  mBeta    <- matrix(0,MAXITER,p)
  vsigma2  <- matrix(0,MAXITER,1)
  vlambda2 <- matrix(0,MAXITER,1)
  #mB <- matrix(0,MAXITER,p)
  #mA <- matrix(0,MAXITER,p)
  
  va_til <- matrix(0,MAXITER,1)
  vb_til <- matrix(0,MAXITER,1)
  vu_til <- matrix(0,MAXITER,1)
  vv_til <- matrix(0,MAXITER,1)
  mM <- matrix(0,MAXITER,p)
  mV <- matrix(0,MAXITER,p)
  
  for (ITER in 1:MAXITER) 
  {
    a_til <- a + 0.5*(n + p)
    b_til <- b + 0.5*sum((vy - mX%*%vbeta)^2) + 0.5*lambda2*sum(vbeta*vbeta/va)
    b_til <- as.vector(b_til)
    #cat("a_til=",a_til,"\n")
    #cat("b_til=",b_til,"\n")
    sigma2 <- invgamma::rinvgamma(1, shape=a_til, rate=b_til)
    
    if ((sigma2<0)|is.na(sigma2)) {
      cat("something has gone horribly wrong with sigma2 \n")
      cat("ITER=",ITER,"\n")
      cat("a_til=",a_til,"\n")
      cat("b_til=",b_til,"\n")
      cat("sum(is.na(vbeta))=",sum(is.na(vbeta)),"\n")
      cat("sum(is.nan(vbeta))=",sum(is.nan(vbeta)),"\n")
      cat("sum(is.na(va))=",sum(is.na(va)),"\n")
      cat("sum(is.nan(va))=",sum(is.nan(va)),"\n")
    }
    sigma <- sqrt(sigma2)
    
    
    if (p<n) {
      mQ_inv <- XTX + diag(lambda2/va)
      mQ <- solve(mQ_inv,tol=1.0E-99)
      mSigma_til <- sigma2*mQ
      vmu_til <- mQ%*%XTy
      vmu_til <- as.vector(vmu_til)
      vbeta <- t(rmvnorm(n = 1, mean = vmu_til, sigma = mSigma_til))
    } else {
      vd  <- lambda2/va
      vd_inv <- 1/vd
      vu <- rnorm(p)*sqrt(vd_inv)
      vv <- mX%*%vu + rnorm(n)
      
      # Upper triangular Cholesky factorization
      # Thhe commented line below is slow
      #mG <- mX%*%diag(vd_inv)%*%t(mX) + mI_n
      
      if (trunc!=0) {
        inds <- which(vd_inv>trunc) 
        mX1 <- mX[,inds]
        mG <- mX1%*%(vd_inv[inds]*t(mX1)) + mI_n
      } else {
        mG <- mX%*%(vd_inv*t(mX)) + mI_n
      }
      
      
      mU <- chol(mG)
      
      mM3 <- backsolve(t(mU), mX, upper.tri=FALSE)
      vv1 <- t(mM3*mM3)%*%matrix(1,n,1) 
      vsigma2_chol <- sigma2*(vd_inv - vd_inv*vd_inv*as.vector(vv1))
      
      vw <- backsolve(mU, backsolve(t(mU), -vv, upper.tri=FALSE))
      vmu <- t(mX)%*%backsolve(mU, backsolve(t(mU), vy, upper.tri = FALSE))
      vmu <- vd_inv*vmu
      
      vbeta <- vmu + sqrt(sigma2)*(vu + vd_inv*(t(mX)%*%vw))
    }
      
    #cat("max(abs(vmu - vmu))",max(abs(vmu_til - vmu)),"\n") 
    #cat("max(abs(diag(mSigma_til) - vsigma2_chol))",max(abs(diag(mSigma_til) - vsigma2_chol)),"\n")  
    #ans <- readline()
    
    vbeta <- as.vector(vbeta)
    

    lambda <- sqrt(lambda2)
    
    #cat("vbeta=",vbeta,"\n")
    
    
    vnu <- sigma/(lambda*abs(vbeta))  
    #vB <- lambda2*vbeta*vbeta/sigma2 
    
 
    va <- 1/rinvgaussian(p,vnu,rep(1,p));
    
    # if (aux_sampler=="statmod") { va = 1/statmod::rinvgauss(p,mean=vnu); } # Can be vectorised
    # if (aux_sampler=="actuar") {  va = 1/actuar::rinvgauss(p,mean=vnu); } # Can be vectorised
    # if (aux_sampler=="LaplacesDemon") { va = 1/LaplacesDemon::rinvgaussian(p,mu=vnu,1); }# Can be vectorised
    # if (aux_sampler=="rmutil_rginvgauss") { va = 1/rmutil::rinvgauss(p,m=vnu, s=1); }# Can be vectorised
    # if (aux_sampler=="rmutil_rinvgauss") { va = rmutil::rginvgauss(p,m=vB, s=1, f=0.5); } # Can be vectorised
    # if (aux_sampler=="SuppDists") { va = 1/SuppDists::rinvGauss(p, nu = vnu, 1); } # Can be vectorised
    # if (aux_sampler=="GIGrvg") { va = GIGrvg::rgig(p,0.5,vB,1); } # Can be vectorised
    # 
    # for (j in 1:p) {
    #   #cat("ITER=",ITER," vB[",j,"]=",vB[j]," sigma2=",sigma2," lambda2=",lambda2," beta[",j,"]=",vbeta[j],"\n")
    #   if (aux_sampler=="GeneralizedHyperbolic") { va[j] = GeneralizedHyperbolic::rgig(1, chi = vB[j], psi = 1, lambda = 0.5); }
    #   if (aux_sampler=="ghyp") { va[j] = ghyp::rgig(1, 0.5, vB[j], 1); }
    # }
    va <- as.vector(va)
    
    if (any(va<0)|any(is.na(va))) {
      cat("something has gone horribly wrong with va \n")
      cat("ITER=",ITER,"\n")
      cat("vnu=",vnu,"\n")
      at("sigma2=",sigma2,"\n")
      at("lambda2=",lambda2,"\n")
      cat("sum(is.na(vbeta))=",sum(is.na(vbeta)),"\n")
      cat("sum(is.nan(vbeta))=",sum(is.nan(vbeta)),"\n")
      cat("sum(is.na(va))=",sum(is.na(va)),"\n")
      cat("sum(is.nan(va))=",sum(is.nan(va)),"\n")
    }
    
    if (USEPRIOR) {
      if (lambda_list$prior=="gamma") {
        u <-   lambda_list$u
        v <-   lambda_list$v
        u_til <- u + 0.5*p
        v_til <- v + 0.5*sum(vbeta*vbeta/va)/sigma2
        lambda2 <- rgamma(1, u_til, v_til)
        
        if ((lambda2<0)|is.na(lambda2)) {
          cat("something has gone horribly wrong with lambda2 \n")
          cat("ITER=",ITER,"\n")
          cat("u_til=",u_til,"\n")
          cat("v_til=",v_til,"\n")
          cat("sum(is.na(vbeta))=",sum(is.na(vbeta)),"\n")
          cat("sum(is.nan(vbeta))=",sum(is.nan(vbeta)),"\n")
          cat("sum(is.na(va))=",sum(is.na(va)),"\n")
          cat("sum(is.nan(va))=",sum(is.nan(va)),"\n")
        }
      } 
    } else {
      u_til <- NA
      v_til <- NA
    }
    
    
    mBeta[ITER,] <- vbeta
    vsigma2[ITER] <- sigma2
    vlambda2[ITER] <- lambda2
    #mB[ITER,] <- vb
    #mA[ITER,] <- va
    
    va_til[ITER] <- a_til
    vb_til[ITER] <- b_til
    vu_til[ITER] <- u_til
    vv_til[ITER] <- v_til
    mM[ITER,] <- vmu_til
    
    if (p<n) {
      mV[ITER,] <- diag(mSigma_til)
    } else {
      mV[ITER,] <- vsigma2_chol
    }
    
    if ((ITER%%1000)==0)
      cat("ITER=",ITER,"\n")
  }
  
  return(
    list(
      mBeta = mBeta,
      vsigma2 = vsigma2,
      vlambda2 = vlambda2,
      va_til = va_til,
      vb_til = vb_til,
      vu_til = vu_til,
      vv_til = vv_til,
      mM = mM,
      mV = mV
    )
  )
}



density_rw_sigma2 <- function(x, va_til, vb_til) 
{
  y <- x*0
  N <- length(va_til)
  for (i in 1:N) {
    y <- y + invgamma::dinvgamma(x, shape=va_til[i], rate=vb_til[i])    
  }
  y <- y/N
  return(list(x=x,y=y))
}

density_rw_lambda2 <- function(x, vu_til, vv_til) 
{
  y <- x*0
  N <- length(vu_til)
  for (i in 1:N) {
    y <- y + dgamma(x, shape=vu_til[i], rate=vv_til[i])  
  }
  y <- y/N
  return(list(x=x,y=y))
}


density_rw_beta <- function(x, mM, mV,j) 
{
  y <- x*0
  N <- nrow(mM)
  for (i in 1:N) {
    y <- y + dnorm(x, mM[i,j], sqrt(mV[i,j]))  
  }
  y <- y/N
  return(list(x=x,y=y))
}





  
  
  