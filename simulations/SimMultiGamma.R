Draw.boxplotsM <- function(Obj)
{
  p0hat <- Obj$p0hat.SP
  Sensitivity <- Obj$Sensitivity
  FPR <- 1 - Obj$Specificity

  par(mfrow = c(1, 3))
  boxplot(p0hat, main = "Estimates of p0")
  boxplot(FPR, main = "False Positive Rate = (1 - Specificity)", ylim = c(0, 0.05))
  boxplot(Sensitivity, main = "Sensitivity", ylim = c(0.8, 1))
}

SimMultGamma <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 2
  
  rho12 <- 0.5
  Sigma0 <- matrix(c(1, rho12, 
                     rho12, 1), 
                   ncol = 2, byrow = T)
  
  library(mvtnorm)
  
  result <- data.frame(p0hat.SP = rep(NA, M),
                       Sensitivity = rep(NA, M),
                       Specificity = rep(NA, M))
  result.1D.1 <- data.frame(p0hat.SP = rep(NA, M),
                            Sensitivity = rep(NA, M),
                            Specificity = rep(NA, M))
  result.1D.2 <- data.frame(p0hat.SP = rep(NA, M),
                            Sensitivity = rep(NA, M),
                            Specificity = rep(NA, M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, frankCopula(param = param, dim = 2))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qgamma(V, shape = 12, rate = 4)
    z <- rbind(z0, z1)
    plot(z, col = c(rep(1, n0), rep(2, n1)), pch = 20)
    
    res <- sp.mix.multi(z)
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    
    res <- sp.mix.1D(z[,1], doplot = FALSE)
    
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result.1D.1$p0hat.SP[r] <- p0hat
    result.1D.1$Sensitivity[r] <- TP/(TP + FN)
    result.1D.1$Specificity[r] <- TN/(TN + FP)
    
    res <- sp.mix.1D(z[,2], doplot = FALSE)
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result.1D.2$p0hat.SP[r] <- p0hat
    result.1D.2$Sensitivity[r] <- TP/(TP + FN)
    result.1D.2$Specificity[r] <- TN/(TN + FP)
  }
  
  return(list(res2D = result, 
              res1D.1 = result.1D.1, 
              res1D.2 = result.1D.2))
}

source(file = '../SpMix.R')

M <- 50
N <- 1000
Res.10 <- SimMultGamma(M = M, n = N, p0 = 0.95)
Res.11 <- SimMultGamma(M = M, n = N, p0 = 0.90)
Res.12 <- SimMultGamma(M = M, n = N, p0 = 0.80)

S <- cbind(Res.10$res1D.1$Sensitivity,
           Res.10$res1D.2$Sensitivity,
           Res.10$res2D$Sensitivity)
colnames(S) <- c("1D", "1D", "2D")
boxplot(S, main = "Sensitivity: p0 = 0.95, f1 = Gamma",
        col = c("skyblue", "skyblue", "orange"), 
        ylim = c(0.5, 1), points = 20)

S <- cbind(Res.11$res1D.1$Sensitivity,
           Res.11$res1D.2$Sensitivity,
           Res.11$res2D$Sensitivity)
colnames(S) <- c("1D", "1D", "2D")
boxplot(S, main = "Sensitivity: p0 = 0.90, f1 = Gamma",
        col = c("skyblue", "skyblue", "orange"), 
        ylim = c(0.5, 1), points = 20)

S <- cbind(Res.12$res1D.1$Sensitivity,
           Res.12$res1D.2$Sensitivity,
           Res.12$res2D$Sensitivity)
colnames(S) <- c("1D", "1D", "2D")
boxplot(S, main = "Sensitivity: p0 = 0.80, f1 = Gamma",
        col = c("skyblue", "skyblue", "orange"), 
        ylim = c(0.5, 1), points = 20)

save.image(file = "SimMultiGamma.RData")

#png(file = "MultiGammap95.png", height = 600, width = 900)
Draw.boxplotsM(Res.10)
#dev.off()

#png(file = "MultiGammap90.png", height = 600, width = 900)
Draw.boxplotsM(Res.11)
#dev.off()

#png(file = "MultiGammap80.png", height = 600, width = 900)
Draw.boxplotsM(Res.12)
#dev.off()
