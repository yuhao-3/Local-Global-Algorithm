###################  MULTIVARIATE LASSO DISTRIBUTION ###############################
library(MASS)
library(mvtnorm)
library(sn)
# Note: A: positive definite,c vb in R2, cc>0, 
# Note: a>0, c>0
dmlasso <- function(x,a,b,c,logarithm=FALSE) 
{
  log_Z <- alasso(a,b,c) 
  log_pdf <-  -0.5*t(x)%*%A%*%x + b %*% x - c*norm(x,"1") - log_Z
  
  if (logarithm) {
    return(log_pdf) 
  }  
  return(exp(log_pdf))
}

# Calculate normalizing constant for Multivariate Lasso 
zmlasso <- function(A,b,c,logarithm=FALSE)
{
  # Get parameter
  A_star = A * matrix(c(1,-1,-1,1),2,2)
  Sigma1 = solve(A)
  Sigma2 = solve(A_star)
  mu1 = as.vector(Sigma1 %*%(b - c))
  mu2 = as.vector(Sigma2 %*% matrix(c(b[1]-c,-b[2]-c),2,1))
  mu3 = as.vector(Sigma2 %*% matrix(c(-b[1]-c,b[2]-c),2,1))
  mu4 = as.vector(Sigma1 %*%(-b - c))
  
  # Get 4 parts of equation separately
  PartI = log(pmvnorm(0,Inf,mean = mu1, sigma = Sigma1))-log(dmvnorm(as.vector(A%*%mu1),sigma = A))
  PartII = log(pmvnorm(0,Inf,mean = mu2, sigma =Sigma2))-log(dmvnorm(as.vector(A_star%*%mu2),sigma = A_star))
  PartIII = log(pmvnorm(0,Inf,mean = mu3, sigma = Sigma2))-log(dmvnorm(as.vector(A_star%*%mu3),sigma = A_star))
  PartIV = log(pmvnorm(0,Inf,mean = mu4, sigma=Sigma1))-log(dmvnorm(as.vector(A%*%mu4),sigma = A))

  # Log Sum Exp Trick to prevent overflow and underflow    
  log_Z1 =  log(det(Sigma2)) + logSumExp(c(PartII[1],PartIII[1]))
  log_Z2 =  log(det(Sigma1)) + logSumExp(c(PartI[1],PartIV[1]))

  if (logarithm) {
    return(log_Z)
  } 
  return(exp(log_Z1) + exp(log_Z2))
}



# Calculate expectation of lasso distribution
emlasso<- function(A,b,c)
{
  Z = zmlasso(A,b,c)
  mean = 
  
  return(mean)
}

# Calculate convariance matrix of multivariate lasso distribution
vmlasso = function(PartI, PartII, PartIII, PartIV, Z)
{
  Sigma =  (PartI*multi_pos_tru_mean(mu1,sigma) + PartII*pos_tru_xxt(mu1,sigma) + PartIII + PartIV)/Z
  return(Sigma)
}

# A = matrix(c(2,0,0,3),2,2)
# b = matrix(c(-2,1),2,1)
# c = 3
# calculate_multi_lasso_param(A,b,c)



# library(pracma)
# fun <- function(x,y) 
# {
#   A = matrix(c(2,0,0,3),2,2)
#   b = matrix(c(-2,1),2,1)
#   c = 3
#   vx = matrix(c(x,y),2,1)
#   return(exp((-0.5*t(vx)%*%A%*%vx + t(b)%*%vx-  c*norm(vx,"1"))[1]))
# }
# xmin <- -10; xmax <- 10
# ymin <- -10; ymax <- 10
# 
# # Normalizing constant
# integral2(fun, xmin, xmax, ymin, ymax)


