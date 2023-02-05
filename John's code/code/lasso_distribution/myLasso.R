###################  UNIVARIATE LASSO DISTRIBUTION ###############################
# Note: a>0,c>0,b in R
calculate_uni_lasso_param <- function(a,b,c)
{
  mu1 = (b-c)/a
  mu2 = -(c+b)/a
  sigma2 = 1/a
  
  Z = norm_const_lasso(mu1,mu2,sqrt(sigma2))
  mu = mean_uni_lasso(mu1,mu2,sqrt(sigma2),Z)
  var = var_uni_lasso(mu1,mu2,sqrt(sigma2),Z)
  return(list(mu = mu,var = var, Z = Z))
}
# Calculate normalizing constant
norm_const_lasso <- function(mu1,mu2,sigma)
{
  # Log trick
  PartI = pnorm(mu1/sigma,log.p = TRUE)- dnorm(mu1/sigma,log = TRUE)
  PartII = pnorm(mu2/sigma,log.p = TRUE)- dnorm(mu2/sigma,log = TRUE)
  
  # Avoid Overflow
  M = max(c(PartI,PartII))
  log_val = M + log(exp(PartI-M) + exp(PartII-M))
  
  # Avoid Underflow
  Z = exp(log(sigma)+log_val)
  return(Z)
}

# Calculate expectation of lasso distribution
mean_uni_lasso <- function(mu1,mu2,sigma,Z)
{
  PartI = (pnorm(mu1/sigma,log.p = TRUE)- dnorm(mu1/sigma,log = TRUE) + log(pos_tru_exp(mu1,sigma)))
  PartII = (pnorm(mu2/sigma,log.p = TRUE)- dnorm(mu2/sigma,log = TRUE) + log(pos_tru_exp(mu2,sigma)))
  
  log_val = PartI - PartII
  
  # Avoid Underflow
  mean = exp(log(sigma) -log(Z) + log_val)
  return(mean)
}


mean_uni_lasso <- function(mu1,mu2,sigma,Z)
{
  part1 = exp(pnorm(mu1/sigma,log.p = TRUE)- dnorm(mu1/sigma,log = TRUE))*pos_tru_exp(mu1,sigma)
  part2 = exp(pnorm(mu2/sigma,log.p = TRUE)- dnorm(mu2/sigma,log = TRUE))*pos_tru_exp(mu2,sigma)
  mean = sigma/Z *(part1 - part2 )
  return(mean)
}

# Calculate second moment of lasso distribution
second_mom_uni_lasso = function(mu1,mu2,sigma,Z)
{
  PartI = (pnorm(mu1/sigma,log.p = TRUE)- dnorm(mu1/sigma,log = TRUE) + log(pos_tru_sec_mom(mu1,sigma)))
  PartII = (pnorm(mu2/sigma,log.p = TRUE)- dnorm(mu2/sigma,log = TRUE) +log(pos_tru_sec_mom(mu2,sigma)))
  
  # Avoid Overflow
  M = max(c(PartI,PartII))  
  log_val = M + log(exp(PartI-M) + exp(PartII-M))
  
  # Avoid Underflow
  sec_mom = exp(log(sigma) -log(Z) + log_val)
  return(sec_mom)
}

# Calculate variance of lasso distribution
var_uni_lasso <- function(mu1,mu2,sigma,Z)
{
  var = second_mom_uni_lasso(mu1,mu2,sigma,Z) -(mean_uni_lasso(mu1,mu2,sigma,Z))^2
  return(var)
}
#################################  Positive Truncated Norm ############################
# Calculate expectation of truncated normal distribution
pos_tru_exp <- function(mu,sigma)
{
  return(mu + sigma*xsi_1(mu/sigma))
}
# Calculate second moment of truncated normal distribution
pos_tru_sec_mom <- function(mu,sigma)
{
  return(sigma^2*(1+xsi_2(mu/sigma)) +(mu + sigma * xsi_1(mu/sigma))^2)
}

# xsi 1 function
xsi_1<- function(t)
{
  return(exp(dnorm(t,log = TRUE)- pnorm(t,log.p = TRUE)))
}
# xsi 2 function
xsi_2<- function(t)
{
  return(-t*xsi_1(t) - xsi_1(t)^2)
}






## John's code

library(sn)
# return the expected value
elasso <- function(a,b,c)
{
  mu1 <-  (b - c)/a
  mu2 <- -(b + c)/a
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  z1 <- zeta(1,mu1/sigma)
  z2 <- zeta(1,mu2/sigma)
  
  e1 <- mu1 + sigma*z1
  e2 <- mu2 + sigma*z2
  
  R1 <- z1/(z1 + z2)
  R2 <- z2/(z1 + z2)
  
  #val <- (e1/z1 - e2/z2)/(1/z1 + 1/z2)
  val <- R2*e1 - R1*e2
  
  return(val)
}
################################################################################
# return the variance
vlasso <- function(a,b,c)
{
  mu1 <-  (b - c)/a
  mu2 <- -(b + c)/a
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  z1 <- zeta(1,mu1/sigma)
  z2 <- zeta(1,mu2/sigma)
  
  e1 <- mu1 + sigma*z1
  e2 <- mu2 + sigma*z2
  
  v1 <- sigma2*(1 + zeta(2,mu1/sigma))
  v2 <- sigma2*(1 + zeta(2,mu2/sigma))
  
  R1 <- z1/(z1 + z2)
  R2 <- z2/(z1 + z2)
  
  #val <- ((v1 + e1*e1)/z1 + (v2 + e2*e2)/z2)/(1/z1 + 1/z2)
  #val <- ((v1 + e1*e1)*z2 + (v2 + e2*e2)*z1)/(z1 + z2)
  val <- R2*(v1 + e1*e1)  + R1*(v2 + e2*e2)
  val <- val - elasso(a,b,c)^2
  
  return(val)
}