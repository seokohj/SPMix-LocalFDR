source("../SpMix.R")

z1 <- read.csv("L.csv")
z2 <- read.csv("O.csv")
z3 <- read.csv("S.csv")

# > dim(z1)
# [1] 536   7
# > dim(z2)
# [1] 415   7
# > dim(z3)
# [1] 503   7

gs <- unique(c(as.character(z1[,1]),
               as.character(z2[,1]),
               as.character(z3[,1])))

id1 <- match(gs, z1[,1])
id2 <- match(gs, z2[,1])
id3 <- match(gs, z3[,1])

id <- cbind(id1, id2, id3)

tmp <- apply(id, 1, is.na)
tmp <- apply(tmp, 2, sum)
idA <- id[!is.na(id[,1]) & !is.na(id[,2])&!is.na(id[,3]),]

z <- cbind(z1[idA[,1],], z2[idA[,2],], z3[idA[,3],])

z0 <- z[,c(5, 12, 19)]
z0 <- -qnorm(as.matrix(z0))
z0[,3] <- pmax(z0[,3], -3)

res <- sp.mix.multi(z0)

res0 <- sp.mix.multi(z0, mono = FALSE) # no monotone constraint

#png(file = "comparison.png", width = 600, height = 600)
plot(res0$localfdr, res$localfdr, 
     xlab = "Local FDR without monotone constraint",
     ylab = "Local FDR with monotone constraint")
abline(0, 1, col = 2, lty = 2)
abline(h = thre, col = 2)
abline(v = thre, col = 2)
#dev.off()


#png(file = "1d_fit.png", width = 900, height = 450)
par(mfrow = c(1, 3))
res1 <- sp.mix.1D(z0[,1])
title("L")
res2 <- sp.mix.1D(z0[,2])
title("O")
res3 <- sp.mix.1D(z0[,3])
title("S")
#dev.off()

thre <- 0.1 # threshold for localFDR
print(BayesFactor <- res$p.0/(1 - res$p.0)*(1/thre - 1)) # Resulting Bayes factor

color <- 1 + as.integer(res$localfdr <= thre) 
color2 <- rep(NA, length(color))
color2[color == 1] <- "gray"
color2[color != 1] <- "black"
pch2 <- rep(1, length(color))
pch2[color != 1] <- 19

library(scatterplot3d)

#png(file = "3dscatter_.png", width = 800, height = 600)
par(mfrow = c(1, 1))
scatterplot3d(x = z0[,1], y = z0[,2], z = z0[,3], 
              color = color2,
              pch = pch2, 
              col.grid = "lightblue",
              xlab = "peripheral leukocytes", 
              ylab = "orbit", 
              zlab = "sinus brushings")
#dev.off()

#library(rgl)
#library(RColorBrewer)
#plot3d(z0[,1], z0[,2], z0[,3], col = color)

z <- cbind(z,
           res$localfdr, 
           res0$localfdr,
           res1$localfdr, 
           res2$localfdr, 
           res3$localfdr,
           res$localfdr <= thre, 
           res0$localfdr <= thre,
           res1$localfdr <= thre, 
           res2$localfdr <= thre, 
           res3$localfdr <= thre)
write.csv(z, file = "localFDR.csv")

save.image("Study1.RData")
