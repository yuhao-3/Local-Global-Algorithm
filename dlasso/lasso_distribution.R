

expit <- function(x) { return( 1/(1+ exp(-x))); }

logSumExp <- function(vec) {
  M <- max(vec)
  val <- M + log(sum(exp(vec - M)))
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


calculate_lasso_dist_stats <- function(a,b,c) 
{
  mu_plus  <-  (b - c)/a
  mu_minus <-  (b + c)/a
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  r_plus  <- mu_plus/sigma
  r_minus <- mu_minus/sigma
  
  z_plus  <- zeta(1, r_plus)
  z_minus <- zeta(1,-r_minus)
  
  log_pm <- pnorm(-(b + c)*sigma, log=TRUE)
  log_pp <- pnorm( (b - c)*sigma, log=TRUE)
  w <- expit(log_pm - log_pp + 2*b*c*sigma2)
  
  
  return(list(
    mu_plus=mu_plus,
    mu_minus=mu_minus,
    sigma2=sigma2,
    sigma=sigma,
    z_plus=z_plus,
    z_minus=z_minus,
    r_plus=r_plus,
    r_minus=r_minus,    
    w=w
  ))	
}


# Note: a>0, c>0
dlasso <- function(x,a,b,c,logarithm=FALSE) 
{
  log_Z <- zlasso(a,b,c,logarithm=TRUE) 
  log_pdf <-  -0.5*a*x*x + b*x - c*abs(x) - log_Z
  
  if (logarithm) {
    return(log_pdf)  
  }  
  return(exp(log_pdf))
}

################################################################################

# return the normalizing constant
zlasso <- function(a,b,c,logarithm=FALSE) 
{
  res <- calculate_lasso_dist_stats(a,b,c)
  log_inv_z_plus  <- pnorm( res$r_plus,log=TRUE)  - dnorm( res$r_plus,log=TRUE)
  log_inv_z_minus <- pnorm(-res$r_minus,log=TRUE) - dnorm(-res$r_minus,log=TRUE)
  log_Z <- log(res$sigma) + logSumExp(c(log_inv_z_plus,log_inv_z_minus))
  if (logarithm) {
    return(log_Z)
  } 
  return(exp(log_Z))
}

################################################################################

plasso <- function(x,a,b,c) 
{
  res <- calculate_lasso_dist_stats(a,b,c)
  
  val <- 0*x
  
  inds_minus <- which(x<=0)
  inds_plus  <- which(x >0)
  
  if (length(inds_minus)>0) {
    val[inds_minus] <- res$w*exp(pnorm((x[inds_minus]-res$mu_minus)/res$sigma,log=TRUE) - pnorm(-res$r_minus,log=TRUE))
  }
  
  if (length(inds_plus)>0) {
    val[inds_plus] <- res$w + (1-res$w)*(1 - exp(pnorm((res$mu_plus - x[inds_plus])/res$sigma, log=TRUE) - pnorm(res$r_plus,log=TRUE)))
    
  }
  
  return(val)
}

################################################################################

# A version of qnorm that doesn't return -Inf when p is very close to 0

qnorm_safe <- function( log_p ) 
{
  TOL <- -500  
  inds1 <- which( log_p <  TOL)
  inds2 <- which( log_p >= TOL)
  
  val <- 0*log_p
  if (length(inds1)>0) {
    val[inds1] <- -sqrt( -2*log_p[inds1] - log(-log_p[inds1]) - log(4*pi) )
  }
  if (length(inds2)>0) {
    val[inds2] <- qnorm(exp(log_p[inds2]))
  }
  return(val)
}

################################################################################


qlasso_fast <- function(u,a,b,c) 
{
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  log_pm <- pnorm(-(b + c)*sigma, log=TRUE)
  log_pp <- pnorm( (b - c)*sigma, log=TRUE)
  
  w <- expit(log_pm - log_pp + 2*b*c*sigma2)
  
  if (u<=w) {
    x <-  (b + c)*sigma2 + sigma*qnorm(exp(log_pm)*u/w)
  } else {
    x <- (b - c)*sigma2 - sigma*qnorm(exp(log_pp)*(1 - u)/(1-w))
  }
  return(x)
}


qlasso_fast4 <- function(u,a,b,c)
  
{
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  r_plus  <- (b - c) * sigma
  r_minus <- (b + c) * sigma
  
  if ( (-4 < r_plus) && (r_plus < 4) && (-4 < r_minus) && (r_minus < 4) ) {
    w <- pnorm(-r_minus) * exp(2 * b * c * sigma2) / pnorm(r_plus)
  } else {
    w <- exp( pnorm(-r_minus, log=TRUE) 
              - pnorm(r_plus, log=TRUE)
              + 2 * b * c * sigma2 )
  }
  
  cat("r_plus=",r_plus,"\n")
  cat("r_minus=",r_minus,"\n")
  cat("w=",w,"\n")
  cat("b=",b,"\n")
  cat("c=",c,"\n")
  cat("sigma2=",sigma2,"\n")
  
  if (u<=w) {
    p1 <- pnorm(-r_minus) * u/w
    mu_minus <-  (b + c) * sigma2
    x <- mu_minus + sigma * qnorm( p1 )
  } else {
    p2 <- pnorm(r_plus) * (1 - u)/(1-w)
    mu_plus  <-  (b - c) * sigma2
    x <- mu_plus - sigma * qnorm( p2 )
  }
  return(x)
}

qlasso <- function(u,a,b,c) 
{
  res <- calculate_lasso_dist_stats(a,b,c)
  
  inds1 <- which(u<=res$w)
  inds2 <- which(u>res$w)
  inds3 <- which(u==0)
  inds4 <- which(u==1)
  
  x <- 0*u
  
  if (length(inds1)>0) {
    log_p1 <- pnorm(-res$r_minus,log=TRUE) + log(u[inds1]) - log(res$w)
    x[inds1]<- res$mu_minus + res$sigma*qnorm_safe( log_p1 )
  }
  
  if (length(inds2)>0) {
    log_p2 <- pnorm(res$r_plus,log=TRUE) + log(1 - u[inds2]) - log(1-res$w)
    x[inds2] <- res$mu_plus - res$sigma*qnorm_safe( log_p2 )
  }
  
  if (length(inds3)>0) {
    x[inds3] <- -Inf
  }
  
  if (length(inds4)>0) {
    x[inds4] <- Inf
  }
  
  
  return(x)
}

################################################################################

rlasso <- function(n,a,b,c) {
  u <- runif(n)
  x <- qlasso(u,a,b,c)
  return(x)
}

rlasso_fast <- function(a,b,c) {
  u <- runif(1)
  x <- qlasso_fast(u,a,b,c)
  return(x)
}

################################################################################

library(sn)

# return the expected value
elasso <- function(a,b,c) 
{
  res <- calculate_lasso_dist_stats(a,b,c)
  z_plus  <- zeta(1, res$r_plus)
  z_minus <- zeta(1,-res$r_minus)
  e_plus  <-  res$mu_plus  + res$sigma*z_plus
  e_minus <-  res$mu_minus - res$sigma*z_minus
  val <- res$w*e_minus + (1-res$w)*e_plus 
  return(val)
}

################################################################################

# return the variance
vlasso <- function(a,b,c) 
{
  res <- calculate_lasso_dist_stats(a,b,c)
  
  z_plus  <- zeta(1, res$r_plus)
  z_minus <- zeta(1,-res$r_minus)
  e_plus  <- res$mu_plus  + res$sigma*z_plus
  e_minus <- res$mu_minus - res$sigma*z_minus
  
  v_plus  <- res$sigma2*(1 + zeta(2,res$r_plus))
  v_minus <- res$sigma2*(1 + zeta(2,-res$r_minus))
  
  val <- res$w*(v_minus + e_minus*e_minus)  + (1-res$w)*(v_plus + e_plus*e_plus)
  val <- val - (res$w*e_minus + (1-res$w)*e_plus)^2
  
  return(val)
}

################################################################################

mlasso <- function(a,b,c) {
  val <- max(c(abs(b)-c,0))*sign(b)/a
  return(val)
}

test_moments <- function(a,b,c) {
  
  # Calculate the moments using the formula
  cat("expected value = ", elasso(a,b,c), "\n")
  cat("variance = ", vlasso(a,b,c) , "\n")
  
  # Calculate moments via trapezoidal integration
  EPS <- 1.0E-6
  L <- qlasso(EPS,a,b,c)
  R <- qlasso(1-EPS,a,b,c)
  N <- 100000
  x <- seq(L,R,length=N)
  f1 <- dlasso(x,a,b,c,logarithm=FALSE) 
  
  zval <- trapint(x,f1)
  eval <- trapint(x,x*f1)
  eval2 <- trapint(x,x*x*f1)
  vval <- eval2 - eval^2
  
  cat("trapint normalizing constant = ", zval, "\n")
  cat("trapint expected value = ", eval, "\n")
  cat("trapint variance = ", vval , "\n")
  
  # Calculate moments via Monte Carlo
  z <- rlasso(n=1000000,a,b,c)
  
  cat("rlasso expected value = ", mean(z), "\n")
  cat("rlasso variance = ", var(z) , "\n")
}