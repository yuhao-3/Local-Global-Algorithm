
################################################################################

# Note: a>0, c>0
dlasso <- function(x,a,b,c,logarithm=FALSE) 
{
	log_Z <- alasso(a,b,c) 
	log_pdf <-  -0.5*a*x*x + b*x - c*abs(x) - log_Z
	
	if (logarithm) {
		return(log_pdf) {
	}  
	return(exp(log_pdf))
}

################################################################################

# return the log-partition function
alasso <- function(a,b,c) 
{	
	Z <- zlasso(a,b,c)
	A <- log(Z)
	return(A)
}

################################################################################

mills_ratio <- function(x) {
	val <- 1/zeta(1,x)
	return(val)
}

################################################################################

# return the normalizing constant
zlasso <- function(a,b,c) 
{
	mu1 <-  (b - c)/a
	mu2 <- -(b + c)/a
	sigma2 <- 1/a
	sigma <- sqrt(sigma2)
	
	Z <- sigma*(mills_ratio(mu1/sigma) + mills_ratio(mu2/sigma))
	return(Z)
}

################################################################################

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
	
	val <- (e1/z1 - e2/z2)/(1/z1 + 1/z2)
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
	
	v1 <- sigma2*(1 + zeta2(2,mu1/sigma))
	v2 <- sigma2*(1 + zeta2(2,mu2/sigma))
	
	val <- ((v1 + e1*e1)/z1 + (v2 + e2*e2)/z2)/(1/z1 + 1/z2)
	val <- val - elasso(a,b,c)^2
	
	return(val)
}


