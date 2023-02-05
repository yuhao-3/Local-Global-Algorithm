dlasso <- function(x,a,b,c) {
  
  mu1 <- (c+b)/a
  mu2 <- (c-b)/a
  sigma2 <- 1/a
  sigma <- sqrt(sigma2)
  
  log_t1 <- 0.5*mu1*mu1/sigma2 + pnorm(mu1/sigma,log=TRUE)
  log_t2 <- 0.5*mu2*mu2/sigma2 + pnorm(mu2/sigma,log=TRUE)
  M <- max(c(log_t1,log_t2))
    
  log_Z <- 0.5*log(2*pi*sigma2) + M + log(exp(log_t1 - M) + exp(log_t2 - M))
  
  #Z <- sqrt(2*pi*sigma2)*(exp((mu1^2)/(2*sigma2))*pnorm(mu1/sigma) + 
  #                          exp((mu2^2)/(2*sigma2))*pnorm(mu2/sigma))
  y <- exp( - 0.5*a*x*x + b*x + c*abs(x) - log_Z) #/Z
  return(y)
}