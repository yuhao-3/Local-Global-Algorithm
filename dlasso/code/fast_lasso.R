################################################################################
# Slightly modified version of rinvgauss function from statmod package
################################################################################

rinvgauss. <- function(n,
                       mean = 1,
                       shape = NULL,
                       dispersion = 1) {
  if (!is.null(shape)) {
    dispersion <- 1 / shape
  }
  mu <- rep_len(mean, n)
  phi <- rep_len(dispersion, n)
  r <- rep_len(0, n)
  i <- (mu > 0 & phi > 0)
  if (!all(i)) {
    r[!i] <- NA
    n <- sum(i)
  }
  phi[i] <- phi[i] * mu[i]
  Y <- rchisq(n, df = 1)
  X1 <- 1 + phi[i]/2*(Y-sqrt(4*Y/phi[i]+Y^2))
  # Note: The line above should yield all X1>0, but it occasionally doesn't due to
  # numerical precision issues. The line below detects this and recomputes
  # the relevant elements of X1 using a 2nd-order Taylor expansion of the
  # sqrt function, which is a good approximation whenever the problem occurs.
  if (any(X1<=0)) {
    X1[X1<=0] <- (1/(Y*phi[i]))[X1<=0]
  }
  firstroot <- as.logical(rbinom(n, size=1L, prob=1/(1+X1)))
  r[i][firstroot] <- X1[firstroot]
  r[i][!firstroot] <- 1/X1[!firstroot]
  mu * r
}

################################################################################
# Draw from inverse-Gaussian distribution while avoiding potential numerical problems
################################################################################

rinvgaussian <- function(n, m, l) {
  m. <- m / sqrt(m * l)
  l. <- l / sqrt(m * l)
  sqrt(m * l) * rinvgauss.(n, m., l.)
}

################################################################################
# Gibbs iteration functions for both Bayesian lassos
# Note: The versions have separate functions, as opposed to being different
# options of the same function, since the latter would require checking any such
# options every time the function is called, i.e., in every MCMC iteration.
# Note: The values of XTY, n, and p can obviously be calculated from X and Y, but they
# are included as inputs to avoid recalculating them every time the function is
# called, i.e., in every MCMC iteration.
################################################################################

iter.bl.original <- function(beta, sigma2, X, Y, XTY, n, p, lambda) 
{
  d.tau.inv <- rinvgaussian(p, sqrt(lambda^2*sigma2/ beta^2), lambda^2)
  d.tau <- as.vector(1/d.tau.inv)
  
  ### Efficient algorithm to generate beta when n < p
  u <- rnorm(p)*sqrt(d.tau)
  v <- X %*% u + rnorm(n)
  U <- chol(X %*% diag(d.tau) %*% t(X) + diag(rep(1, n)))
  w <- backsolve(U, backsolve(t(U), ((Y/sqrt(sigma2)) - v), upper.tri=F))
  
  beta.new <- sqrt(sigma2)*(u + d.tau*(t(X) %*% w))
  sigma2.new <- (sum((Y-drop(X%*%beta.new))^2) + sum(beta.new^2*d.tau.inv))/rchisq(1,n+p-1)
  
  return(list(beta = beta.new, sigma2 = sigma2.new))
}

################################################################################

iter.bl.fast <- function(beta, sigma2, X, Y, XTY, n, p, lambda) 
{
  d.tau.inv <- rinvgaussian(p, sqrt(lambda^2*sigma2/beta^2), lambda^2)
  d.tau <- as.vector(1/d.tau.inv)
  
  u <- rnorm(p)*sqrt(d.tau)
  v <- X%*% u + rnorm(n)
  U <- chol(X %*%diag(d.tau)%*%t(X) + diag(rep(1, n)))
  w <- backsolve(U, backsolve(t(U), -v, upper.tri = F))
  
  ### Efficient algorithm to generate sigma when n < p
  beta.tilde <- t(X)%*%backsolve(U, backsolve(t(U), Y, upper.tri = F))
  beta.tilde <- d.tau*beta.tilde
  sigma2.new <- (sum(Y^2) - sum(XTY*beta.tilde))/rchisq(1, n - 1)
  
  ### Efficient algorithm to generate beta when n < p
  beta.new <- beta.tilde + sqrt(sigma2.new)*(u + d.tau*(t(X)%*%w))
  
  return(list(beta = beta.new, sigma2 = sigma2.new))
}

################################################################################
# Run original and fast Bayesian lassos
################################################################################

run.bl <-
  function(X,
           Y,
           lambda,
           K,     # number of iterations
           M,     # Number of chains
           outfile.stem,
           fast = F,
           keep.beta = T,
           write.each = F) {
    XTY <- drop(t(X) %*% Y)
    n <- dim(X)[1]
    print(n)
    p <- dim(X)[2]
    print(p)
    iter.bl <- get(paste("iter.bl.", ifelse(fast, "fast", "original"), sep = ""))
    chaindigits <- max(1, ceiling(log(M, 10))) # digits needed for chain label strings
    for (chain in 0:(M - 1)) {
      beta <- rep(1, p)
      sigma2 <- 1
      chaintext <- substring(format(chain / (10 ^ chaindigits), nsmall = chaindigits), 3)
      outfile.beta <- paste(outfile.stem, "-", chaintext, "-b.txt", sep = "")
      outfile.sigma2 <- paste(outfile.stem, "-", chaintext, "-s.txt", sep = "")
      if (write.each) {
        for (k in 1:K) {
          iter.result <- iter.bl(beta, sigma2, X, Y, XTY, n, p, lambda)
          beta <- iter.result$beta
          sigma2 <- iter.result$sigma2
          if(keep.beta) {
            beta.row <- matrix(beta, nrow = 1)
            write.table(
              beta.row,
              outfile.beta,
              append = T,
              row.names = F,
              col.names = F
            )
          }
          write.table(
            sigma2,
            outfile.sigma2,
            append = T,
            row.names = F,
            col.names = F
          )
        }
      } else{
        beta.chain <- matrix(NA, nrow = K, ncol = p)
        sigma2.chain <- rep(NA, K)
        for (k in 1:K) {
          iter.result <- iter.bl(beta, sigma2, X, Y, XTY, n, p, lambda)
          beta <- iter.result$beta
          sigma2 <- iter.result$sigma2
          beta.chain[k, ] <- beta
          sigma2.chain[k] <- sigma2
          
          if ((k%%10)==0) {
            cat("chain=",chain,"k=",k,"\n")
          }
        }
        if (keep.beta) {
          write.table(beta.chain,
                      outfile.beta,
                      row.names = F,
                      col.names = F)
        }
        write.table(sigma2.chain,
                    outfile.sigma2,
                    row.names = F,
                    col.names = F)
      }
      typetext <- ifelse(fast, "Fast", "Orig")
      print(paste(typetext, "chain", chain + 1, "of", M, "complete at", date()))
      flush.console()
    }
    
    return(list(beta.chain = beta.chain, sigma2.chain=sigma2.chain))
  }

################################################################################

# Gibbs iteration functions for both spike-slab samplers
# Note: The versions have separate functions, as opposed to being different
# options of the same function, since the latter would require checking any such
# options every time the function is called, i.e., in every MCMC iteration.
# Note: The values of XTY, n, and p can obviously be calculated from X and Y, but they
# are included as inputs to avoid recalculating them every time the function is
# called, i.e., in every MCMC iteration.

iter.ss.original <- function(beta, sigma2, X, Y, XTY, n, p, w, kappa, zeta) 
{
  w.t <- 1/(1+(1/w-1)*sqrt(kappa)*exp(beta^2*(1/kappa-1)/(2*sigma2*zeta)))
  d.tau <- zeta * (1 + (kappa - 1) * rbinom(p, 1, w.t))
  
  ### Efficient algorithm to generate beta when n < p
  u <- rnorm(p) * sqrt(d.tau)
  v <- X %*% u + rnorm(n)
  U <- chol(X %*% diag(d.tau) %*% t(X) + diag(rep(1, n)))
  ww <- backsolve(U, backsolve(t(U), ((Y/sqrt(sigma2))-v), upper.tri=F))
  beta.new <- sqrt(sigma2) * (u + d.tau*(t(X)%*%ww))
  sigma2.new <- (sum((Y-drop(X%*%beta.new))^2) + sum(beta.new^2/d.tau))/rchisq(1, n+p-1)
  
  return(list(beta = beta.new, sigma2 = sigma2.new))
}

################################################################################

iter.ss.fast <- function(beta, sigma2, X, Y, XTY, n, p, w, kappa, zeta) 
{
  w.t <- 1/(1+(1/w-1)*sqrt(kappa)*exp(beta^2*(1/kappa-1)/(2*sigma2*zeta)))
  d.tau <- zeta * (1 + (kappa - 1) * rbinom(p, 1, w.t))
  u <- rnorm(p) * sqrt(d.tau)
  v <- X %*% u + rnorm(n)
  U <- chol(X %*% diag(d.tau) %*% t(X) + diag(rep(1, n)))
  ww <- backsolve(U, backsolve(t(U), -v, upper.tri = F))
  
  ### Efficient algorithm to generate sigma when n < p
  beta.tilde <- t(X) %*% backsolve(U, backsolve(t(U), Y, upper.tri = F))
  beta.tilde <- d.tau * beta.tilde
  sigma2.new <- (sum(Y ^ 2) - sum(XTY * beta.tilde)) / rchisq(1, n - 1)
  
  ### Efficient algorithm to generate beta when n < p
  beta.new <- beta.tilde + sqrt(sigma2.new)*(u + d.tau*(t(X)%*%ww))
  
  return(list(beta = beta.new, sigma2 = sigma2.new))
}

################################################################################
# Run original and fast spike-slab samplers
################################################################################

run.ss <-
  function(X,
           Y,
           w,
           kappa,
           zeta,
           K,
           M,
           outfile.stem,
           fast = F,
           keep.beta = T,
           write.each = F) {
    XTY <- drop(t(X) %*% Y)
    n <- dim(X)[1]
    p <- dim(X)[2]
    iter.ss <-
      get(paste("iter.ss.", ifelse(fast, "fast", "original"), sep = ""))
    chaindigits <-
      max(1, ceiling(log(M, 10))) # digits needed for chain label strings
    for (chain in 0:(M - 1)) {
      beta <- rep(1, p)
      sigma2 <- 1
      chaintext <-
        substring(format(chain / (10 ^ chaindigits), nsmall = chaindigits), 3)
      outfile.beta <-
        paste(outfile.stem, "-", chaintext, "-b.txt", sep = "")
      outfile.sigma2 <-
        paste(outfile.stem, "-", chaintext, "-s.txt", sep = "")
      if (write.each) {
        for (k in 1:K) {
          iter.result <- iter.ss(beta, sigma2, X, Y, XTY, n, p, w, kappa, zeta)
          beta <- iter.result$beta
          sigma2 <- iter.result$sigma2
          if (keep.beta) {
            beta.row <- matrix(beta, nrow = 1)
            write.table(
              beta.row,
              outfile.beta,
              append = T,
              row.names = F,
              col.names = F
            )
          }
          write.table(
            sigma2,
            outfile.sigma2,
            append = T,
            row.names = F,
            col.names = F
          )
        }
      } else{
        beta.chain <- matrix(NA, nrow = K, ncol = p)
        sigma2.chain <- rep(NA, K)
        for (k in 1:K) {
          iter.result <- iter.ss(beta, sigma2, X, Y, XTY, n, p, w, kappa, zeta)
          beta <- iter.result$beta
          sigma2 <- iter.result$sigma2
          beta.chain[k, ] <- beta
          sigma2.chain[k] <- sigma2
        }
        if (keep.beta) {
          write.table(beta.chain,
                      outfile.beta,
                      row.names = F,
                      col.names = F)
        }
        write.table(sigma2.chain,
                    outfile.sigma2,
                    row.names = F,
                    col.names = F)
      }
      typetext <- ifelse(fast, "Fast", "Orig")
      print(paste(typetext, "chain", chain + 1, "of", M, "complete at", date()))
      flush.console()
    }
  }        

################################################################################