Draw.boxplots <- function(Obj)
{
  Sensitivity.SP <- Obj$TP.SP/(Obj$TP.SP + Obj$FN.SP)
  Sensitivity.EF <- Obj$TP.EF/(Obj$TP.EF + Obj$FN.EF)
  Specificity.SP <- Obj$TN.SP/(Obj$TN.SP + Obj$FP.SP)
  Specificity.EF <- Obj$TN.EF/(Obj$TN.EF + Obj$FP.EF)
  
  p0hat <- with(Obj, cbind(p0hat.EF, p0hat.SP))
  RMSE <- with(Obj, cbind(RMSE.EF, RMSE.SP))
  Sensitivity <- cbind(Sensitivity.EF, Sensitivity.SP)
  FPR <- cbind(1-Specificity.EF, 1-Specificity.SP)
  colnames(p0hat) <- colnames(RMSE) <- colnames(Sensitivity) <- colnames(FPR) <- c("Efron", "Proposed")
  
  par(mfrow = c(2, 2))
  boxplot(p0hat, main = "Estimates of p0")
  boxplot(RMSE, main = "RMSE of local FDR estimates")
  boxplot(FPR, main = "False Positive Rate = (1 - Specificity)", ylim = c(0, 0.05))
  boxplot(Sensitivity, main = "Sensitivity", ylim = c(0, 1))
}

SimModel2 <- function(M, n, p0, alpha, beta)
  # [Univariate] Null: Normal, Nonnull: Gamma
{
  library(locfdr)
  
  result <- data.frame(p0hat.SP = rep(NA, M),
                       RMSE.SP = rep(NA, M),
                       TP.SP = rep(NA, M),
                       TN.SP = rep(NA, M),
                       FP.SP = rep(NA, M),
                       FN.SP = rep(NA, M),
                       p0hat.EF = rep(NA, M),
                       RMSE.EF = rep(NA, M),
                       TP.EF = rep(NA, M),
                       TN.EF = rep(NA, M),
                       FP.EF = rep(NA, M),
                       FN.EF = rep(NA, M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    z0 <- rnorm(n0)
    z1 <- rgamma(n1, shape = alpha, rate = (1/beta))
    z <- c(z0, z1)
    tmp <- p0*dnorm(z)
    f <- tmp + (1-p0)*dgamma(z, shape = alpha, rate = (1/beta))
    localFDR <- tmp/f
    
    res <- sp.mix.1D(z, doplot = TRUE)
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    RMSE <- sqrt(mean((localFDR >= 0.01)*(localFDR <= 0.5)*(res$localfdr - localFDR)^2))
    result$p0hat.SP[r] <- p0hat
    result$RMSE.SP[r] <- RMSE
    result$TP.SP[r] <- TP
    result$TN.SP[r] <- TN
    result$FP.SP[r] <- FP
    result$FN.SP[r] <- FN
    
    res <- locfdr(z)
    p0hat <- res$fp0[3, 3]
    Nhat <- as.integer(res$fdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    RMSE <- sqrt(mean((localFDR >= 0.01)*(localFDR <= 0.5)*(res$fdr - localFDR)^2))
    result$p0hat.EF[r] <- p0hat
    result$RMSE.EF[r] <- RMSE
    result$TP.EF[r] <- TP
    result$TN.EF[r] <- TN
    result$FP.EF[r] <- FP
    result$FN.EF[r] <- FN
  }
  
  return(result)
}

#source(file = '../SpMix.R')

#Res.4 <- SimModel2(M = 500, p0 = 0.95, n = 1000, alpha = 12, beta = .25)
#Res.5 <- SimModel2(M = 500, p0 = 0.90, n = 1000, alpha = 12, beta = .25)
#Res.6 <- SimModel2(M = 500, p0 = 0.80, n = 1000, alpha = 12, beta = .25)
#save.image(file = "SimUniGamma.RData")

load(file = "SimUniGamma.RData")

#png(file = "UniGammap95.png", height = 600, width = 900)
Draw.boxplots(Res.4)
#dev.off()

#png(file = "UniGammap90.png", height = 600, width = 900)
Draw.boxplots(Res.5)
#dev.off()

#png(file = "UniGammap80.png", height = 600, width = 900)
Draw.boxplots(Res.6)
#dev.off()

