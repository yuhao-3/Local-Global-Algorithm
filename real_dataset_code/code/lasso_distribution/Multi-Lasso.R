###################  MULTIVARIATE LASSO DISTRIBUTION ###############################
library(MASS)
library(mvtnorm)
library(sn)
library(matrixStats)
# library(TruncatedNormal)
# library(hpa)
# library(truncnorm)
library(mnormt)

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
    return(log(exp(log_Z1) + exp(log_Z2)))
  } 
  return(exp(log_Z1) + exp(log_Z2))
}



# Calculate expectation of multivariate lasso distribution
emlasso<- function(A,b,c)
{
  # Get parameter
  lower = 0
  upper = Inf
  A_star = A * matrix(c(1,-1,-1,1),2,2)
  Sigma1 = solve(A)
  Sigma2 = solve(A_star)
  mu1 = as.vector(Sigma1 %*%(b - c))
  mu2 = as.vector(Sigma2 %*% matrix(c(b[1]-c,-b[2]-c),2,1))
  mu3 = as.vector(Sigma2 %*% matrix(c(-b[1]-c,b[2]-c),2,1))
  mu4 = as.vector(Sigma1 %*%(-b - c))
  vsigma1 = Sigma1[row(Sigma1)==col(Sigma1)]
  vsigma2 = Sigma2[row(Sigma2)==col(Sigma2)]
  
  
  # Get normalizign constant
  Z = zmlasso(A,b,c)
  # Calculate expectation for each of the four parts
  tn_mean1 = log(mom.mtruncnorm(powers = 1, mu1, Sigma1, rep(lower,2), rep(upper,2))$cum1)
  tn_mean2 = log(mom.mtruncnorm(powers = 1, mu2, Sigma2, rep(lower,2), rep(upper,2))$cum1)
  tn_mean3 = log(mom.mtruncnorm(powers = 1, mu3, Sigma2, rep(lower,2), rep(upper,2))$cum1)
  tn_mean4 = log(mom.mtruncnorm(powers = 1, mu4, Sigma1, rep(lower,2), rep(upper,2))$cum1)

  # # Get 4 parts of equation separately
  # PartI = log(pmvnorm(0,Inf,mean = mu1, sigma = Sigma1)) + tn_mean1-log(dmvnorm(as.vector(A%*%mu1),sigma = A))
  # PartII =(log(pmvnorm(0,Inf,mean = mu2, sigma = Sigma2)) + tn_mean2-log(dmvnorm(as.vector(A_star%*%mu2),sigma = A_star)))*c(1,-1) 
  # PartIII = (log(pmvnorm(0,Inf,mean = mu3, sigma = Sigma2)) + tn_mean3-log(dmvnorm(as.vector(A_star%*%mu3),sigma = A_star)))*c(-1,1)
  # PartIV = (log(pmvnorm(0,Inf,mean = mu4, sigma = Sigma1)) + tn_mean4-log(dmvnorm(as.vector(A%*%mu4),sigma = A))) *c(-1,-1)
  # 
  # 
  # 
  # # Log Sum Exp Trick to prevent overflow and underflow
  # log_mean1 =  log(det(Sigma2)) + log(exp(PartII)+ exp(PartIII)) - log(Z)
  # log_mean2 =  log(det(Sigma1)) + log(exp(PartI)+ exp(PartIV)) - log(Z)
  # 
  # mean = exp(log_mean1)+ exp(log_mean2)
  # 
  # # 
  # # 
  mean = det(Sigma1)/Z * (exp(tn_mean1)* pmvnorm(0,Inf,mean = mu1, sigma = Sigma1)
                          /dmvnorm(as.vector(A%*%mu1),sigma = A)
                        +c(-1,-1)*exp(tn_mean4)* pmvnorm(0,Inf,mean = mu4, sigma = Sigma1)
                        /dmvnorm(as.vector(A%*%mu4),sigma = A)) +
         det(Sigma2)/Z *(c(1,-1)*exp(tn_mean2)* pmvnorm(0,Inf,mean = mu2, sigma = Sigma2)
                                  /dmvnorm(as.vector(A_star%*%mu2),sigma = A_star)
                        +c(-1,1)*exp(tn_mean3)* pmvnorm(0,Inf,mean = mu3, sigma = Sigma2)
                        /dmvnorm(as.vector(A_star%*%mu3),sigma = A_star))
  
  


  
  return(mean)
}

# Calculate convariance matrix of multivariate lasso distribution
vmlasso = function(A,b,c)
{
  # Get parameter
  A_star = A * matrix(c(1,-1,-1,1),2,2)
  Sigma1 = solve(A)
  Sigma2 = solve(A_star)
  mu1 = as.vector(Sigma1 %*%(b - c))
  mu2 = as.vector(Sigma2 %*% matrix(c(b[1]-c,-b[2]-c),2,1))
  mu3 = as.vector(Sigma2 %*% matrix(c(-b[1]-c,b[2]-c),2,1))
  mu4 = as.vector(Sigma1 %*%(-b - c))
  
  # # Get 4 parts of equation separately
  # PartI = log(pmvnorm(0,Inf,mean = mu1, sigma = Sigma1))-log(dmvnorm(as.vector(A%*%mu1),sigma = A))
  # PartII = log(pmvnorm(0,Inf,mean = mu2, sigma =Sigma2))-log(dmvnorm(as.vector(A_star%*%mu2),sigma = A_star))
  # PartIII = log(pmvnorm(0,Inf,mean = mu3, sigma = Sigma2))-log(dmvnorm(as.vector(A_star%*%mu3),sigma = A_star))
  # PartIV = log(pmvnorm(0,Inf,mean = mu4, sigma=Sigma1))-log(dmvnorm(as.vector(A%*%mu4),sigma = A))
  # 
  
  Z = zmlasso(A,b,c)
  mu = emlasso(A,b,c)

  
  # Calculate expectation for each of the four parts
  tn_cov1 = mom.mtruncnorm(powers = 2, mu1, Sigma1, rep(lower,2), rep(upper,2))$order2$m2
  tn_cov2 = mom.mtruncnorm(powers = 2, mu2, Sigma2, rep(lower,2), rep(upper,2))$order2$m2
  tn_cov3 = mom.mtruncnorm(powers = 2, mu3, Sigma2, rep(lower,2), rep(upper,2))$order2$m2
  tn_cov4 = mom.mtruncnorm(powers = 2, mu4, Sigma1, rep(lower,2), rep(upper,2))$order2$m2
  
  
  
  EX2 = det(Sigma1)/Z * (tn_cov1* pmvnorm(0,Inf,mean = mu1, sigma = Sigma1)
                          /dmvnorm(as.vector(A%*%mu1),sigma = A)
                          +tn_cov4* pmvnorm(0,Inf,mean = mu4, sigma = Sigma1)
                          /dmvnorm(as.vector(A%*%mu4),sigma = A)) +
        det(Sigma2)/Z *(matrix(c(1,-1,-1,1),2,2)*tn_cov2* pmvnorm(0,Inf,mean = mu2, sigma = Sigma2)
                    /dmvnorm(as.vector(A_star%*%mu2),sigma = A_star)
                    +matrix(c(1,-1,-1,1),2,2)*tn_cov3* pmvnorm(0,Inf,mean = mu3, sigma = Sigma2)
                    /dmvnorm(as.vector(A_star%*%mu3),sigma = A_star))

  print(EX2)

  cov = EX2 - mu %*% t(mu)
  
  return(cov)
  
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


