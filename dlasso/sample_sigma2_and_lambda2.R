################################################################################

sample_sigma2 <- function(sigma2, a_val, b_val, c_val) 
{
  MAXITER <- 1000
  TOL <- 1.0E-6
  PLOT <- FALSE
  inflation_factor <- 1.5
  
  m <- sigma2
  for (ITER in 1:MAXITER) {
    #f_val <- -a_val*log(m) - b_val/m - c_val/sqrt(m)
    g_val <- -a_val/m + b_val/(m*m) + 0.5*c_val/(sqrt(m)*m)
    h_val <- a_val/(m*m) - 2*b_val/(m*m*m) - 1.5*0.5*c_val/(sqrt(m)*m*m)
    m <- m - g_val/h_val
    
    if (m<0) {
      m <- 1.0E-12
    }
    
    #cat(ITER,m,   f_val, g_val, h_val, "\n")
    
    if (abs(g_val)<TOL) {
      break;
    }
  }
  
  log_p_fun <- function(sigma2, a_val, b_val, c_val) {
    val <- -a_val*log(sigma2) - b_val/sigma2 - c_val/sqrt(sigma2)
  }
  
  s2 <- -1/h_val
  
  s2 <- s2*inflation_factor
  
  a_til <- m*m/s2 + 2
  b_til <- m*(a_til - 1)
  
  if (PLOT) {
    s <- sqrt(s2)
    nsd <- 5
    N <- 100
    x <- seq(max(c(0,m-nsd*s)),m+nsd*s,length=N)
    f <- log_p_fun(x, a_val, b_val, c_val)
    f <- exp(f - max(f))
    f <- f/trapint(x,f)
    plot(x,f,type="l")
    lines(x,dinvgamma(x,shape=a_til, rate=b_til),col="red")
    ans <- readline()
  }
  
  sigma2_prop <- rinvgamma(1, shape=a_til, rate=b_til)
  
  log_p_prop <- log_p_fun(sigma2_prop, a_val, b_val, c_val)
  log_q_prop <- dinvgamma(sigma2_prop, shape=a_til, rate=b_til,log=TRUE)
  log_p_curr <- log_p_fun(sigma2, a_val, b_val, c_val)
  log_q_curr <- dinvgamma(sigma2, shape=a_til, rate=b_til, log=TRUE)
  
  alpha <- exp(log_p_prop - log_p_curr +  log_q_curr - log_q_prop)
  alpha <- min(c(1,alpha))
  
  accept <- FALSE
  if (runif(1) < alpha) {
    sigma2 <- sigma2_prop
    accept <- TRUE
  }
  
  return(list(sigma2=sigma2,accept=accept))
}

################################################################################

slice_sampler <- function(x, a_val, b_val, c_val) 
{
  PLOT <- FALSE
  TOL <- 1.0E-8
  MAXITER <- 100
  
  # The target distribution
  f_fun <- function(x, a_val, b_val, c_val) {
    val <- a_val*log(x) - b_val*x - c_val*sqrt(x);
    return(val);
  }
  
  # It's graadient
  g_fun <- function(x, a_val, b_val, c_val) {
    val <- a_val/x - b_val - 0.5*c_val/sqrt(x);
    return(val);
  }
  
  # It's Hessian
  h_fun <- function(x, a_val, b_val, c_val) {
    val <- -a_val/(x*x)  + 0.25*c_val/(sqrt(x)*x)
    return(val);
  }
  
  # Calculate the mode
  z_star <- (-0.5*c_val + sqrt(0.25*c_val*c_val + 4*a_val*b_val))/(2*b_val)
  x_star <- z_star*z_star
  
  # Target at the mode
  f_star <- f_fun(x_star, a_val, b_val, c_val)
  p_star <- exp(f_star)
  
  # Calculate approximate variaiance and sd
  h_star <- h_fun(x_star, a_val, b_val, c_val) 
  sigma2 <- -1/h_star
  sigma  <- sqrt(sigma2)
  
  # Calculate target at initial x
  x0 <- x
  f0 <- f_fun(x0, a_val, b_val, c_val) - f_star
  p0 <- exp(f0)
  
  if (p0==0) {
    # Something has gone horribly wrong because the 
    # initial x is far outside the effective domain
    x0 <- x_star
    f0 <- f_fun(x0, a_val, b_val, c_val) - f_star
    p0 <- exp(f0)
  }
  
  # Sample uniform vertically
  u <- runif(1,0,p0)
  log_u <- log(u)

  # Initial guesses of the left and right end-points based on Laplace approximation
  xL <- x_star - sqrt(2*sigma2*( - log_u))
  xR <- x_star + sqrt(2*sigma2*( - log_u))

  if (PLOT) {
    cat("x=", x, "\n")
    cat("a_val=", a_val, "\n")
    cat("b_val=", b_val, "\n")
    cat("c_val=", c_val, "\n")
    cat("z_star=", z_star, "\n")
    cat("x_star=", x_star, "\n")
    cat("f_star=", f_star, "\n")
    cat("p_star=", p_star, "\n")
    cat("h_star=", h_star, "\n")
    cat("sigma2=", sigma2, "\n")
    cat("sigma=",  sigma, "\n")
    cat("x0=", x0, "\n")
    
    cat("f0=", f0, "\n")
    cat("p0=", p0, "\n")
    cat("u=", u, "\n")
    cat("log_u=", log_u, "\n")
    cat("xL=", xL, "\n")
    cat("xR=", xR, "\n")
  }
    
  if (PLOT) {
    vx <- seq(xL-sigma,xR+sigma,length=300)
    vf <- f_fun(vx, a_val, b_val, c_val) - f_star
    vp <- exp(vf)

    plot(vx, vp, type="l")
    points(x_star, f_star)
    lines(rep(x_star,2), c(0, p_star), col="red")
    vn <- p_star*exp( -((vx - x_star)^2)/(2*sigma2))
    lines(vx,vn, col="blue")
    lines(rep(x0,2), c(0, p0), col="red")
    points(x0, u)
    
    lines(c(xL, xR), rep(u,2),col="green")
  }
    
  # Set left end-point initial guess to be small if Laplace approximation 
  # gives negative end-point
  if (xL<0) {
    xL <- 1.0E-12
  }
  
  # Use Newton's method to find left end point
  for (ITER in 1:MAXITER) {
    fL <- f_fun(xL, a_val, b_val, c_val) - f_star
    gL <- g_fun(xL, a_val, b_val, c_val)
    xL <- xL - (fL - log_u)/gL
    if (abs(gL)<TOL) { break; }
  }
  
  # Use Newton's method to find right end-point
  for (ITER in 1:MAXITER) {
    fR <- f_fun(xR, a_val, b_val, c_val) - f_star
    gR <- g_fun(xR, a_val, b_val, c_val)
    xR <- xR - (fR - log_u)/gR
    if (abs(gR)<TOL) { break; }
  }
  
  # Sample uniformly between left and right end-points
  x_new <- runif(1, xL, xR)
  
  if (PLOT) {
    lines(c(xL, xR), rep(u,2),col="green")
  }
    
  return(x_new)
}


################################################################################

sample_lambda2 <- function(lambda2, a_val, b_val, c_val) 
{
  MAXITER <- 1000
  TOL <- 1.0E-6
  PLOT <- FALSE
  inflation_factor <- 1.5
  
  m <- lambda2
  for (ITER in 1:MAXITER) {
    #f_val <- a_val*log(m) - b_val*m - sqrt(m)*c_val
    g_val <- a_val/m - b_val - 0.5*c_val/sqrt(m)
    h_val <- -a_val/(m*m)  + 0.25*c_val/(sqrt(m)*m)
    m <- m - g_val/h_val
    
    if (m<0) {
      m <- 1.0E-12
    }
    
    #cat(ITER,m,   f_val, g_val, h_val, "\n")
    
    if (abs(g_val)<TOL) {
      break;
    }
  }
  
  log_p_fun <- function(lambda2, a_val, b_val, c_val) {
    val <- a_val*log(lambda2) - b_val*lambda2 - sqrt(lambda2)*c_val
  }
  
  s2 <- -1/h_val
  
  s2 <- s2*inflation_factor
  
  u_til <- m*m/s2
  v_til <- m/s2
  
  if (PLOT) {
    s <- sqrt(s2)
    nsd <- 5
    N <- 100
    x <- seq(max(c(0,m-nsd*s)),m+nsd*s,length=N)
    f <- log_p_fun(x, a_val, b_val, c_val)
    f <- exp(f - max(f))
    f <- f/trapint(x,f)
    plot(x,f,type="l")
    lines(x,dgamma(x,shape=u_til, rate=v_til),col="red")
  }
  
  lambda2_prop <- rgamma(1, shape=u_til, rate=v_til)
  
  log_p_prop <- log_p_fun(lambda2_prop, a_val, b_val, c_val)
  log_q_prop <- dgamma(lambda2_prop, shape=u_til, rate=v_til,log=TRUE)
  log_p_curr <- log_p_fun(lambda2, a_val, b_val, c_val)
  log_q_curr <- dgamma(lambda2, shape=u_til, rate=v_til, log=TRUE)
  
  alpha <- exp(log_p_prop - log_p_curr +  log_q_curr - log_q_prop)
  alpha <- min(c(1,alpha))
  
  accept <- FALSE
  if (runif(1) < alpha) {
    lambda2 <- lambda2_prop
    accept <- TRUE
  }
  
  return(list(lambda2=lambda2,accept=accept))
}
