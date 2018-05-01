sp.mix.logconc <- function(z, doplot = TRUE, tol = 1.0e-4, max.iter = 100)
# For carina data only
{
  library(LogConcDEAD)

  ## Initial step
  q0 <- quantile(z, probs=.7)
  p.0 <- mean(z >= q0)
  mu.0 <- mean(z[z >= q0])
  sig.0 <- sd(z[z >= q0])
  
  mu.1 <- mean(z[z < q0])
  sig.1 <- sd(z[z < q0])
  f1.tilde <- dnorm(z, mean=mu.1, sd=sig.1)
  
  term1 <- p.0*dnorm(z, mu.0, sig.0)
  term2 <- term1 + (1-p.0)*f1.tilde
  gam <- term1/term2
  
  ## EM-step  
  diff <- 100
  k <- 0
  while ((k < 3) | ((k < max.iter)&(diff > tol))) {  
    k <- k+1
    ## E-step
    term1 <- p.0*dnorm(z, mu.0, sig.0)
    term2 <- term1 + (1-p.0)*f1.tilde
    new.gam <- term1/term2
    
    ## M-step
    w.gam <- new.gam/sum(new.gam, na.rm=T)
    new.mu.0 <- sum(w.gam*z, na.rm=T)
    new.sig.0 <- sqrt(sum(w.gam*(z-new.mu.0)^2, na.rm=T))
    new.p.0 <- mean(new.gam, na.rm=T)
    
    which.z <- new.gam <= .9
    weight <- (1-new.gam[which.z])/sum(1-new.gam[which.z])
    f1.tilde[which.z] <- exp(mlelcd(z[which.z], w=weight)$logMLE)
    f1.tilde[!which.z] <- 0
    which.gam <- (new.gam<=0.5)*(new.gam>=0.01)
    diff <- max(abs(gam-new.gam)[which.gam])
    
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0; gam <- new.gam
  }
  
  if (doplot) {
    which.z <- (gam <= .2) & (z <= 250)
    nc <- 100
    thre <- max(z[which.z])
    
    hist(z, probability=T, nclass=nc, xlab="", main="",
        sub = substitute(
            paste(p[0], " = ", p0, ", ", 
                  mu[0], " = ", mu0, ", ",
                  sigma[0], " = ", sigma0,
                  sep = ""), 
            list(p0 = round(p.0, 2),
                 mu0 = round(mu.0, digits = 2),
                 sigma0 = round(sig.0, digits = 2))))
    rug(z, col="gray")
    rug(z[which.z], col=2)
    zs <- sort(z)
    lines(zs, p.0*dnorm(zs, mean=mu.0, sd=sig.0), col=3, lwd = 2)
    lines(zs, (1-p.0)*f1.tilde[order(z)], col=2, lwd = 2)
  }
  res <- list(p.0=p.0, mu.0=mu.0, tau.0=sig.0, localfdr=gam, thre=thre, iter=k)
  
  return(res)
}

dat <- read.table('carina.dat')
str(dat)

x <- dat[dat$V8 + dat$V9 > 0,]
x <- x[x$V6 < 3,]
vel <- x$V4 # represents the Radial velocity data of stars in the Carina galaxy

vel.fit <- sp.mix.logconc(vel)
