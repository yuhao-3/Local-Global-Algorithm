###################  MULTIVARIATE LASSO DISTRIBUTION ###############################
library(MASS)
library(mvtnorm)
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
zmlasso <- function(A,b,c)
{
  # Get parameter
  A_inv = solve(A)
  mu1 = as.vector(A_inv %*%(b - c))
  mu2 = as.vector(A_inv %*% matrix(c(b[1]-c,-b[2]-c),2,1))
  mu3 = as.vector(A_inv %*% matrix(c(-b[1]-c,b[2]-c),2,1))
  mu4 = as.vector(A_inv %*%(-b - c))
  Sigma = A_inv
  
  # Get 4 parts of equation separately
  PartI = pmvnorm(0,Inf,mean = mu1, sigma = Sigma)/dmvnorm(as.vector(A%*%mu1),sigma = A)
  PartII = pmvnorm(0,Inf,mean = mu2, sigma =Sigma)/dmvnorm(as.vector(A%*%mu2),sigma = A)
  PartIII = pmvnorm(0,Inf,mean = mu3, sigma = Sigma)/dmvnorm(as.vector(A%*%mu3),sigma = A)
  PartIV = pmvnorm(0,Inf,mean = mu4, sigma=Sigma)/dmvnorm(as.vector(A%*%mu4),sigma = A)
    
  print(PartIV[1]*det(Sigma))
  
  return((PartI+PartII+PartIII+PartIV)[1]*det(Sigma))
}

# Calculate expectation of lasso distribution
emlasso<- function(A,b,c)
{
  mean = (PartI*multi_pos_tru_mean(mu1,sigma) + PartII*pos_tru_xxt(mu1,sigma) + PartIII + PartIV)/Z
  
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


