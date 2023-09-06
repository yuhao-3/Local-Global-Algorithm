################################################################################

library(statmod)
N <- 100
gauss_quad <- gauss.quad(N, kind = "legendre")

################################################################################

expit <- function(x) { return( 1/(1+exp(-x))); }

################################################################################

logSumExp <- function(vec) {
  M <- max(vec)
  val <- M + log(sum(exp(vec-M)))
  return(val)
}

################################################################################

trapint <- function(xgrid, fgrid) 
{
  ng <- length(xgrid)
  xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
  fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
  integ <- sum(xvec * fvec)/2
  return(integ)
}

################################################################################

library(statmod)

log_pbvnorm_quadrant <- function(vz,rho,quadrant="mm") 
{
  if (quadrant=="pp") {
    return(log_pbvnorm(-vz, rho))
  }
  if (quadrant=="mm") {
    return(log_pbvnorm(vz, rho))
  }
  if (quadrant=="pm") {
    return(log_pbvnorm(c(-vz[1],vz[2]), -rho))
  }
  if (quadrant=="mp") {
    return(log_pbvnorm(c(vz[1],-vz[2]), -rho))
  }
}

pbvnorm_quadrant <- function(vz,mSigma,quadrant="mm") 
{
  vz <- as.vector(vz)
  if (quadrant=="pp") {
    return(pmvnorm(mean=vz, sigma=mSigma, lower=0, upper=Inf))
  }
  if (quadrant=="mm") {
    return(pmvnorm(mean=vz, sigma=mSigma, lower=-Inf, upper=0))
  }
  if (quadrant=="pm") {
    return(pmvnorm(mean=vz, sigma=mSigma, lower=c(0,-Inf), upper=c(Inf,0)))
  }
  if (quadrant=="mp") {
    return(pmvnorm(mean=vz, sigma=mSigma, lower=c(-Inf,0), upper=c(0,Inf)))
  }
}


################################################################################

# Function to calculate the bivariate normal cdf on the log-scale

log_pbvnorm <- function(vz,rho)
{
  N <- 10001
  #cat("vz=",vz," rho=",rho,"\n")
  
  # Obtain the quadrature weights
  #vw <-gauss_quad$weights
  
  # Scale the quadrature points so that they are between 0 and rho
  vt <- seq(0,rho,length=N)  #*rho*sign(rho) runif(N)*rho*sign(rho) 0.5*rho*(gauss_quad$nodes + 1)
  #vt <- 0.5*rho*(gauss_quad$nodes + 1)
  
  vdelta <- abs(diff(vt)[1])
  vw <- rep(2,N)
  vw[1] <- 1
  vw[N] <- 1
  
  # Calculate Gauss-Legendre terms on log-scale using Owen's (1956) integral 
  # formula for the bivariate normal cdf
  #v0 <- sum(pnorm(vz,log=TRUE))
  vv <-   -log(2*pi) - 0.5*log(1 - vt*vt) - 0.5*(vz[1]*vz[1] + vz[2]*vz[2] - 2*vz[1]*vz[2]*vt)/(1 - vt*vt)
  #vv <- c(v0,vv)
  
  #plot(vt,exp(vv),type="l")

  #print(vv)
  #print(log(vw))
  #print(vdelta)
  log_int <- logSumExp(vv + log(vw) + log(vdelta))
  
  # Apply the log-sum-exp function
  val <- logSumExp(c(sum(pnorm(vz,log=TRUE)), log_int)) #+  trapint(vt,exp(vv))
  #val <- log(val)
  #val <- logSumExp(vv)
  
  
  return(val)
}

################################################################################

# Function to calculate various statistics for the bivariate lasso distribution

sbvlasso <- function(mA_val,vb_val,c_val) 
{
  # Calculate the common scale matrix
  mSigma <- solve(mA_val, tol=1.0E-99)  
  
  sigma1 <- sqrt(mSigma[1,1])
  sigma2 <- sqrt(mSigma[2,2])
  vsigma <- c(sigma1,sigma2)
  rho <- mSigma[1,2]/(sigma1*sigma2)
  
  vmu_pp <- mSigma%*%(vb_val - c_val)
  vmu_pm <- mSigma%*%(vb_val - c_val*c(1,-1))
  vmu_mp <- mSigma%*%(vb_val - c_val*c(-1,1))
  vmu_mm <- mSigma%*%(vb_val + c_val)
  
  log_det_mA <- determinant(mA_val)$modulus
  attributes(log_det_mA) <- NULL
  const <-  -0.5*log_det_mA + log(2*pi)
  
  eta_pp <- log(pbvnorm_quadrant(vmu_pp,mSigma,quadrant="pp")) + 0.5*t(vmu_pp)%*%mA_val%*%vmu_pp + const
  eta_pm <- log(pbvnorm_quadrant(vmu_pm,mSigma,quadrant="pm")) + 0.5*t(vmu_pm)%*%mA_val%*%vmu_pm + const
  eta_mp <- log(pbvnorm_quadrant(vmu_mp,mSigma,quadrant="mp")) + 0.5*t(vmu_mp)%*%mA_val%*%vmu_mp + const
  eta_mm <- log(pbvnorm_quadrant(vmu_mm,mSigma,quadrant="mm")) + 0.5*t(vmu_mm)%*%mA_val%*%vmu_mm + const
  veta <- c(eta_pp,eta_pm,eta_mp,eta_mm)
  
  veta_til <- veta - max(veta) 
  vw <- exp(veta_til)/sum(exp(veta_til))
  
  return(list(
    mSigma = mSigma,
    vsigma = vsigma,
    rho = rho,
    vmu_pp = vmu_pp,
    vmu_pm = vmu_pm,
    vmu_mp = vmu_mp,
    vmu_mm = vmu_mm,
    veta = veta,
    vw = vw
  ))
}

################################################################################

# Function to calculate the normalizing constant of the bivariate lasso distribution

zbvlasso <- function(mA_val,vb_val,c_val, logarithm=FALSE) 
{
  res <- sbvlasso(mA_val,vb_val,c_val)
  log_Z <- logSumExp(res$veta) 
  if (logarithm) {
    return(log_Z)
  } 
  return(exp(log_Z))
}

################################################################################

# Calculate the expecation of the bivariate lasso

library(tmvtnorm)

ebvlasso <- function(mA_val,vb_val,c_val,verbose=FALSE) 
{
  res <- sbvlasso(mA_val,vb_val,c_val)
  vw <- res$vw
  
  #print(res)
  
  # These 4 lines could be faster if we used the result specific to the 
  # bivariate truncated normal distribution 
  ve_pp <- mtmvnorm(as.vector(res$vmu_pp), res$mSigma, doComputeVariance=FALSE, 
                    lower=c(0,0), upper=c(Inf,Inf))$tmean
  ve_pm <- mtmvnorm(as.vector(res$vmu_pm), res$mSigma, doComputeVariance=FALSE, 
                    lower=c(0,-Inf), upper=c(Inf,0))$tmean
  ve_mp <- mtmvnorm(as.vector(res$vmu_mp), res$mSigma, doComputeVariance=FALSE, 
                    lower=c(-Inf,0), upper=c(0,Inf))$tmean
  ve_mm <- mtmvnorm(as.vector(res$vmu_mm), res$mSigma, doComputeVariance=FALSE, 
                    lower=c(-Inf,-Inf), upper=c(0,0))$tmean
  
  if (verbose) {
    print("ebvlasso")
    cat("vw=",vw,"\n")
    cat("vmu_pp=",res$vmu_pp,"\n")
    cat("vmu_pm=",res$vmu_pm,"\n")
    cat("vmu_mp=",res$vmu_mp,"\n")
    cat("vmu_mm=",res$vmu_mm,"\n")
    cat("mSigma=",res$mSigma,"\n")
    
    cat("ve_pp=",ve_pp,"\n")
    cat("ve_pm=",ve_pm,"\n")
    cat("ve_mp=",ve_mp,"\n")
    cat("ve_mm=",ve_mm,"\n")
  }  
  
  TOL <- 1.0E-6
  val <- rep(0,length(ve_pp))
  if (vw[1]>TOL) { #a
    val <- val + vw[1]*ve_pp
  }
  if (vw[2]>TOL) { #b
    val <- val + vw[2]*ve_pm
  }
  if (vw[3]>TOL) { #c
    val <- val + vw[3]*ve_mp
  }
  if (vw[4]>TOL) { #d
    val <- val + vw[4]*ve_mm
  }
  
  if (verbose) {
    cat("val=",val,"\n")
  } 

  
  return(val)  
}


################################################################################

# Calculate the variance of the bivariate lasso

vbvlasso <- function(mA_val,vb_val,c_val,verbose=FALSE) 
{
  res <- sbvlasso(mA_val,vb_val,c_val)
  vw <- res$vw
  
  d <- length(res$vmu_pp)
  
  #print(res)
  
  # These 4 lines could be faster if we used the result specific to the 
  # bivariate truncated normal distribution 
  res_pp <- mtmvnorm(as.vector(res$vmu_pp), res$mSigma, doComputeVariance=TRUE, 
                     lower=c(0,0), upper=c(Inf,Inf))
  res_pm <- mtmvnorm(as.vector(res$vmu_pm), res$mSigma, doComputeVariance=TRUE, 
                     lower=c(0,-Inf), upper=c(Inf,0))
  res_mp <- mtmvnorm(as.vector(res$vmu_mp), res$mSigma, doComputeVariance=TRUE, 
                     lower=c(-Inf,0), upper=c(0,Inf))
  res_mm <- mtmvnorm(as.vector(res$vmu_mm), res$mSigma, doComputeVariance=TRUE, 
                     lower=c(-Inf,-Inf), upper=c(0,0))
  
  ve_pp <- res_pp$tmean
  ve_pm <- res_pm$tmean
  ve_mp <- res_mp$tmean
  ve_mm <- res_mm$tmean
  
  ve2_pp <- res_pp$tvar + ve_pp%*%t(ve_pp)
  ve2_pm <- res_pm$tvar + ve_pm%*%t(ve_pm)
  ve2_mp <- res_mp$tvar + ve_mp%*%t(ve_mp)
  ve2_mm <- res_mm$tvar + ve_mm%*%t(ve_mm)
  
  
  if (verbose) {
    print("vbvlasso")
    cat("vw=",vw,"\n")
    cat("vmu_pp=",res$vmu_pp,"\n")
    cat("vmu_pm=",res$vmu_pm,"\n")
    cat("vmu_mp=",res$vmu_mp,"\n")
    cat("vmu_mm=",res$vmu_mm,"\n")
    cat("mSigma=",res$mSigma,"\n")
    
    cat("ve_pp=",ve_pp,"\n")
    cat("ve_pm=",ve_pm,"\n")
    cat("ve_mp=",ve_mp,"\n")
    cat("ve_mm=",ve_mm,"\n")
    cat("ve2_pp=",ve2_pp,"\n")
    cat("ve2_pm=",ve2_pm,"\n")
    cat("ve2_mp=",ve2_mp,"\n")
    cat("ve2_mm=",ve2_mm,"\n")
  }
  
  TOL <- 1.0E-6
  ve <- rep(0,d)
  if (vw[1]>TOL) { #e
    ve <- ve + vw[1]*ve_pp
  }
  if (vw[2]>TOL) { #f
    ve <- ve + vw[2]*ve_pm
  }
  if (vw[3]>TOL) { #g
    ve <- ve + vw[3]*ve_mp
  }
  if (vw[4]>TOL) { #h
    ve <- ve + vw[4]*ve_mm
  }
  ve <- matrix(ve,d,1)
  
  ve2 <- rep(0,d,d)
  if (vw[1]>TOL) ve2 <- ve2 + vw[1]*ve2_pp
  if (vw[2]>TOL) ve2 <- ve2 + vw[2]*ve2_pm
  if (vw[3]>TOL) ve2 <- ve2 + vw[3]*ve2_mp
  if (vw[4]>TOL) ve2 <- ve2 + vw[4]*ve2_mm
  ve2 <- matrix(ve2,d,d)
 
  mV <- ve2 - ve%*%t(ve)
  
  if (verbose) {
    cat("ve=",ve,"\n")
    cat("ve2=",ve2,"\n")
    cat("mV=",mV,"\n")
  }

  
  return(mV)  
}


dmmlasso1 <- function(x,A,b,c,logarithm = FALSE)
{
  sigma2 = 1/A[2,2]
  Z <- zbvlasso(A,b,c)
  k = exp(-0.5*A[1,1]*x^2+b[1]*x-c*abs(x))/Z
  mu1 = -0.5*((A[1,2]+A[2,1])/A[2,2])*x + (b[2]-c)/A[2,2]
  mu2 = 0.5*((A[1,2]+A[2,1])/A[2,2])*x - (b[2]+c)/A[2,2]
  
  log_partI  <- pnorm(mu1/sqrt(sigma2),log=TRUE)  - dnorm(mu1/sqrt(sigma2),log=TRUE)
  log_partII <- pnorm(mu2/sqrt(sigma2),log=TRUE) - dnorm(mu2/sqrt(sigma2),log=TRUE)
  
  M <- max(c(log_partI,log_partII))
  log_mpdf =  log(sqrt(sigma2)) + log(k) + M +log(exp(log_partI - M)+exp(log_partII - M))
  
  if (logarithm) {
    return(log_mpdf)
  }
  return(exp(log_mpdf))
  
}


# Get marginal distribution of x_2
dmmlasso2 <- function(x,A,b,c,logarithm = FALSE)
{
  sigma2 = 1/A[1,1]
  Z <- zbvlasso(A,b,c)
  k = exp(-0.5*A[2,2]*x^2+b[2]*x-c*abs(x))/Z
  mu1 = -0.5*((A[1,2]+A[2,1])/A[1,1])*x + (b[1]-c)/A[1,1]
  mu2 = 0.5*((A[1,2]+A[2,1])/A[1,1])*x - (b[1]+c)/A[1,1]
  
  log_partI  <- pnorm(mu1/sqrt(sigma2),log=TRUE)  - dnorm(mu1/sqrt(sigma2),log=TRUE)
  log_partII <- pnorm(mu2/sqrt(sigma2),log=TRUE) - dnorm(mu2/sqrt(sigma2),log=TRUE)
  
  M <- max(c(log_partI,log_partII))
  log_mpdf =  log(sqrt(sigma2)) + log(k) + M + log(exp(log_partI - M)+exp(log_partII - M))
  
  if (logarithm) {
    return(log_mpdf)
  }
  return(exp(log_mpdf))
  
}
