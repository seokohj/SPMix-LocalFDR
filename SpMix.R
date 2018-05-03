normal.mixture <- function(z, tol = 5e-3, max.iter = 50)
{
  library(mvtnorm)
  
  k <- 0; diff <- 100
  z <- as.matrix(z)
  if ( dim(z)[2] == 1 ) m.dist <- z/sd(z) else m.dist <- mahalanobis(z, center = rep(0, dim(z)[2]), cov = cov(z))
  
  p.0 <- mean(m.dist <= 1.65)
  mu.0 <- rep(0, dim(z)[2])
  sig.0 <- diag(1, dim(z)[2])
  f.0 <- dmvnorm(z, mean = mu.0, sigma = sig.0)
  
  mu.1 <- apply(z[m.dist > 1.65,], 2, mean)
  sig.1 <- cov(z[m.dist > 1.65,])
  f.1 <- dmvnorm(z, mean = mu.1, sigma = sig.1)
  
  while ( (k < 3)|((k < max.iter) & (diff > tol)) ) {
    k <- k + 1
    
    ## E-step
    term1 <- p.0*f.0
    term2 <- term1 + (1-p.0)*f.1
    gam <- term1/term2
    
    ## M-step
    new.p.0 <- mean(gam)
    new.mu.0 <- as.vector(t(z)%*%gam)/sum(gam)
    dev <- (z-new.mu.0)*sqrt(gam)
    new.sig.0 <- t(dev)%*%dev/sum(gam)
    f.0 <- dmvnorm(z, mean = new.mu.0, sigma = new.sig.0)
    new.mu.1 <- as.vector(t(z)%*%(1-gam))/sum(1-gam)
    dev <- (z-new.mu.1)*sqrt(1-gam)
    new.sig.1 <- t(dev)%*%dev/sum(1-gam)
    f.1 <- dmvnorm(z, mean = new.mu.1, sigma = new.sig.1)
    
    ## Update
    diff <- max(abs((new.mu.0 - mu.0)), 
                abs((new.sig.0 - sig.0)), 
                abs((new.mu.1 - mu.1)), 
                abs((new.sig.1 - sig.1)), 
                abs((new.p.0 - p.0)))
    p.0 <- new.p.0
    mu.0 <- new.mu.0
    sig.0 <- new.sig.0
    mu.1 <- new.mu.1
    sig.1 <- new.sig.1
  }
  return(list(iter = k, p.0 = p.0, mu.0 = mu.0, Sigma.0 = sig.0, mu.1 = mu.1, Sigma.1 = sig.1))
}


sp.mix.multi <- function(z, tol = 5e-4, max.iter = 30, mono = TRUE)
  # FOR MULTIVARIATE CASE ONLY
{
  library(LogConcDEAD)
  
  z <- as.matrix(z)
  n <- dim(z)[1]
  
  ## Initial step: to fit normal mixture
  nmEM <- normal.mixture(z)
  p.0 <- nmEM$p.0
  mu.0 <- nmEM$mu.0
  sig.0 <- nmEM$Sigma.0
  f1.tilde <- dmvnorm(z, nmEM$mu.1, nmEM$Sigma.1)
  gam <- f <- rep(0, n)
  
  ## EM-step    
  k <- 0; converged <- 0
  while ( (k < 3)|((k < max.iter) & (!converged)) ) {   
    k <- k + 1
    ## E-step
    tmp <- p.0*dmvnorm(z, mu.0, sig.0)
    new.f <- tmp + (1-p.0)*f1.tilde
    new.gam <- tmp/new.f
    if(mono) new.gam <- MonotoneFDR(z, new.gam)
    
    ## M-step
    sum.gam <- sum(new.gam)
    new.mu.0 <- as.vector(t(z)%*%new.gam)/sum.gam
    dev <- t(t(z)-new.mu.0)*sqrt(new.gam)
    new.sig.0 <- t(dev)%*%dev/sum.gam
    new.p.0 <- mean(new.gam)
    new.f.0 <- dmvnorm(z, new.mu.0, new.sig.0)
    weight <- 1 - new.gam
    new.f1.tilde <- rep(0, n)
    which.z <- (new.gam <= .9)
    lcd <- mlelcd(z[which.z,], w = weight[which.z]/sum(weight[which.z]))
    new.f1.tilde[which.z] <- exp(lcd$logMLE)
    
    ## Update
    which.gam <- (new.gam <= 0.9)*(new.gam >= 0.01)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in mdfdr fit = ", round(diff, 5), "\n") 
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde
    gam <- new.gam 
    f <- new.f
  }
  
  res <- list(p.0 = p.0, mu.0 = mu.0, tau.0 = sig.0, 
              f1.hat = f1.tilde, f = f, localfdr = gam, iter = k)
  
  return(res)
}

sp.mix.1D <- function(z, tol = 5.0e-6, max.iter = 10, doplot = TRUE, thre.localFDR = 0.2)
{
  library(LogConcDEAD)
  #library(fmlogcondens)
  
  z <- as.numeric(z)
  n <- length(z)
  
  ## Initial step
  q0 <- quantile(z, probs = .9)
  p.0 <- mean(z <= q0)
  mu.0 <- mean(z[z <= q0])
  sig.0 <- sd(z[z <= q0])
  
  mu.1 <- mean(z[z > q0])
  sig.1 <- sd(z[z > q0])
  f1.tilde <- dnorm(z, mean = mu.1, sd = sig.1)
  
  f <- gam <- rep(0, n)
  
  ## EM-step  
  k <- 0; converged <- 0
  while ( (k < 3) | ((k < max.iter) & (!converged)) ) {   
    k <- k + 1
    ## E-step
    tmp <- p.0*dnorm(z, mu.0, sig.0)
    new.f <- tmp + (1 - p.0)*f1.tilde
    new.gam <- tmp/new.f
    
    ## M-step
    w.gam <- new.gam/sum(new.gam, na.rm = TRUE)
    new.mu.0 <- sum(w.gam*z, na.rm = TRUE)
    new.sig.0 <- sqrt(sum(w.gam*(z-new.mu.0)^2, na.rm = TRUE))
    new.p.0 <- mean(new.gam, na.rm = TRUE)
    
    new.f1.tilde <- rep(0, n)
    which.z <- new.gam <= .95
    weight <- 1 - new.gam[which.z]
    weight <- weight/sum(weight)
    new.f1.tilde[which.z] <- exp(mlelcd(z[which.z], w = weight)$logMLE)
    #new.f1.tilde[which.z] <- exp(fmlcd(matrix(z[which.z], 
    #                                          nrow = length(weight), 
    #                                          ncol = 1), 
    #                                   w = weight[which.z]/sum(weight[which.z]))$logMLE)
    
    which.gam <- (new.gam <= 0.9)*(new.gam >= 0.01)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in 1dfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde; gam <- new.gam; f <- new.f
  }
  
  which.z <- gam <= thre.localFDR
  thre <- min(z[which.z])
  
  if (doplot) {
    hist(z,
         nclass = max(round(length(z)/20), 24),
         probability = TRUE,
         col = "gray", border = "white", 
         xlab = "", 
         main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ", 
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0, ", ",
                 "threshold = ", threshold, 
                 sep = ""), 
           list(p0 = round(p.0, 2),
                mu0 = round(mu.0, digits = 2),
                sigma0 = round(sig.0, digits = 2),
                threshold = round(thre, digits = 2))))
    rug(z, col = "gray")
    rug(z[which.z], col = 2)
    zs <- sort(z)
    lines(zs, p.0*dnorm(zs, mean = mu.0, sd = sig.0), col = 3, lwd = 2)
    lines(zs, (1-p.0)*f1.tilde[order(z)], col = 2, lwd = 2)
    points(thre, 0, bg = "yellow", col = 2, pch = 25)
  }
  
  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0, f = f,
              localfdr = gam, iter = k)
  
  return(res)
}

NE <- function(x, X)
{
  n <- nrow(X)
  p <- ncol(X)
  xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
  ne.ind <- apply(1*(X >= xx), 1, prod) 
  
  return((1:n)[ne.ind == 1])
}


MonotoneFDR <- function(z, fdr)
{
  n <- nrow(z)
  MFDR <- numeric(n)
  for (i in 1:n) {
    MFDR[i] <- max(fdr[NE(z[i,], z)])
  }
  
  return(MFDR)
}
