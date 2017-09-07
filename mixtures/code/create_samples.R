## full
rm(list = ls())

res <- list()
for (cc in 1:10) {
    res[[cc]] <- readRDS(paste0("/Shared/ssrivastva/wasp/mixtures/result/full/res_", cc, "_100k.rds"))
}

ordMat <- matrix(NA, 10, 2)
for (cc in 1:10) {
    if (all(rowMeans(res[[cc]]$mu[1, , ]) < rowMeans(res[[cc]]$mu[2, , ]))) {
        ordMat[cc, ] <- c(1, 2)
    } else {
        ordMat[cc, ] <- c(2, 1)
    }
}

fullCorr <- list()
for (cc in 1:10) {
    fullCorr[[cc]] <- list(numeric(1000), numeric(1000))
    for (ss in 1:1000) {
        if (all(ordMat[cc, ] == 1:2)) {
            fullCorr[[cc]][[1]][ss] <- cov2cor(res[[cc]]$cov[1, , , ss])[1, 2]
            fullCorr[[cc]][[2]][ss] <- cov2cor(res[[cc]]$cov[2, , , ss])[1, 2]
        } else {
            fullCorr[[cc]][[2]][ss] <- cov2cor(res[[cc]]$cov[1, , , ss])[1, 2]
            fullCorr[[cc]][[1]][ss] <- cov2cor(res[[cc]]$cov[2, , , ss])[1, 2]
        }
    }
}

muList <- list()
for (cc in 1:10) {
    muList[[cc]] <- list()
    if (all(ordMat[cc, ] == 1:2)) {
        muList[[cc]][[1]] <- t(res[[cc]]$mu[1, , ])
        muList[[cc]][[2]] <- t(res[[cc]]$mu[2, , ])
    } else {
        muList[[cc]][[2]] <- t(res[[cc]]$mu[1, , ])
        muList[[cc]][[1]] <- t(res[[cc]]$mu[2, , ])
    }
}

saveRDS(muList, "~/wasp/mixtures/result/muList.rds")
saveRDS(fullCorr, "~/wasp/mixtures/result/corrList.rds")

library(mvtnorm)
xx <- seq(0, 10, length = 500)
densSampMCMC <- rep(list(matrix(NA, 1000, 500)), 10)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    for (ii in 1:1000) {
        cat("ii: ", ii, "\n")
        for (gg in seq_along(xx)) {
            yy <- c(xx[gg], xx[gg])
            dens1 <- dmvnorm(yy, mean = res[[cc]]$mu[1, , ii], sigma = res[[cc]]$cov[1, , , ii])
            dens2 <- dmvnorm(yy, mean = res[[cc]]$mu[2, , ii], sigma = res[[cc]]$cov[2, , , ii])
            densSampMCMC[[cc]][ii, gg] <- res[[cc]]$prob[ii, 1] * dens1 + res[[cc]]$prob[ii, 2] * dens2
        }
    }
}

saveRDS(densSampMCMC, "/Shared/ssrivastva/wasp/mixtures/result/mcmcDensList.rds")

## wasp

rm(list=ls())
setwd("/Shared/ssrivastva/wasp/mixtures/result/sub5/samp")

nsub <- 5
res <- list()
for (cc in 1:10) {
    res[[cc]] <- list()
    for (kk in 1:nsub) {
        res[[cc]][[kk]] <- readRDS(paste0("res_cv_", cc,"_nsub_", kk, "_k5.rds"))
    }
}

muList <- list()
rhoList <- list()
for (cc in 1:10) {
    muList[[cc]] <- list()
    rhoList[[cc]] <- list()
    for (kk in 1:nsub) {
        muList[[cc]][[kk]] <- vector("list", 2)
        names(muList[[cc]][[kk]]) <- c("(1, 2)", "(7, 8)")
        rhoList[[cc]][[kk]] <- vector("list", 2)
        names(rhoList[[cc]][[kk]]) <- c("(1, 2)", "(7, 8)")
    }
}

for (cc in 1:10) {
    for (kk in 1:nsub) {
        if (all(rowMeans(res[[cc]][[kk]]$mu[1, , ]) < rowMeans(res[[cc]][[kk]]$mu[2, , ]))) {
            muList[[cc]][[kk]][[1]] <- t(res[[cc]][[kk]]$mu[1, , ])
            muList[[cc]][[kk]][[2]] <- t(res[[cc]][[kk]]$mu[2, , ])
            rhoList[[cc]][[kk]][[1]] <- numeric(1000)
            rhoList[[cc]][[kk]][[2]] <- numeric(1000)
            for (ss in 1:1000) {
                rhoList[[cc]][[kk]][[1]][ss] <- cov2cor(res[[cc]][[kk]]$cov[1, , , ss])[1, 2]
                rhoList[[cc]][[kk]][[2]][ss] <- cov2cor(res[[cc]][[kk]]$cov[2, , , ss])[1, 2]
            }
        } else {
            muList[[cc]][[kk]][[1]] <- t(res[[cc]][[kk]]$mu[2, , ])
            muList[[cc]][[kk]][[2]] <- t(res[[cc]][[kk]]$mu[1, , ])
            rhoList[[cc]][[kk]][[1]] <- numeric(1000)
            rhoList[[cc]][[kk]][[2]] <- numeric(1000)
            for (ss in 1:1000) {
                rhoList[[cc]][[kk]][[1]][ss] <- cov2cor(res[[cc]][[kk]]$cov[2, , , ss])[1, 2]
                rhoList[[cc]][[kk]][[2]][ss] <- cov2cor(res[[cc]][[kk]]$cov[1, , , ss])[1, 2]
            }
        }
    }
}

saveRDS(rhoList, "~/wasp/mixtures/result/subRhoList5.rds")
saveRDS(muList, "~/wasp/mixtures/result/subMuList5.rds")

for (cc in 1:10) {
    for (kk in 1:nsub) {
        write.table(muList[[cc]][[kk]][[1]], file = paste0("csv/samp_cv_", cc, "_nsub_", kk, "_k5_mu1.csv"),
                    sep = ",", row.names = FALSE, col.names = FALSE)
        write.table(muList[[cc]][[kk]][[2]], file = paste0("csv/samp_cv_", cc, "_nsub_", kk, "_k5_mu2.csv"),
                    sep = ",", row.names = FALSE, col.names = FALSE)
    }
}

library(mvtnorm)
xx <- seq(0, 10, length = 500)
densSampK5 <- rep(list(rep(list(matrix(NA, 1000, 500)), 5)), 10)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    for (kk in 1:nsub) {
        for (ii in 1:1000) {
            cat("ii: ", ii, "\n")
            for (gg in seq_along(xx)) {
                yy <- c(xx[gg], xx[gg])
                dens1 <- dmvnorm(yy, mean = res[[cc]][[kk]]$mu[1, , ii], sigma = res[[cc]][[kk]]$cov[1, , , ii])
                dens2 <- dmvnorm(yy, mean = res[[cc]][[kk]]$mu[2, , ii], sigma = res[[cc]][[kk]]$cov[2, , , ii])
                densSampK5[[cc]][[kk]][ii, gg] <- res[[cc]][[kk]]$prob[ii, 1] * dens1 + res[[cc]][[kk]]$prob[ii, 2] * dens2
            }
        }
    }
}

saveRDS(densSampK5, "/Shared/ssrivastva/wasp/mixtures/result/waspDensList_k5.rds")

rm(list=ls())

setwd("/Shared/ssrivastva/wasp/mixtures/result/sub10/samp")
nsub <- 10

res <- list()
for (cc in 1:10) {
    res[[cc]] <- list()
    for (kk in 1:nsub) {
        res[[cc]][[kk]] <- readRDS(paste0("res_cv_", cc,"_nsub_", kk, "_k10.rds"))
    }
}

muList <- list()
rhoList <- list()
for (cc in 1:10) {
    muList[[cc]] <- list()
    rhoList[[cc]] <- list()
    for (kk in 1:nsub) {
        muList[[cc]][[kk]] <- vector("list", 2)
        names(muList[[cc]][[kk]]) <- c("(1, 2)", "(7, 8)")
        rhoList[[cc]][[kk]] <- vector("list", 2)
        names(rhoList[[cc]][[kk]]) <- c("(1, 2)", "(7, 8)")
    }
}

for (cc in 1:10) {
    for (kk in 1:nsub) {
        if (all(rowMeans(res[[cc]][[kk]]$mu[1, , ]) < rowMeans(res[[cc]][[kk]]$mu[2, , ]))) {
            muList[[cc]][[kk]][[1]] <- t(res[[cc]][[kk]]$mu[1, , ])
            muList[[cc]][[kk]][[2]] <- t(res[[cc]][[kk]]$mu[2, , ])
            rhoList[[cc]][[kk]][[1]] <- numeric(1000)
            rhoList[[cc]][[kk]][[2]] <- numeric(1000)
            for (ss in 1:1000) {
                rhoList[[cc]][[kk]][[1]][ss] <- cov2cor(res[[cc]][[kk]]$cov[1, , , ss])[1, 2]
                rhoList[[cc]][[kk]][[2]][ss] <- cov2cor(res[[cc]][[kk]]$cov[2, , , ss])[1, 2]
            }
        } else {
            muList[[cc]][[kk]][[1]] <- t(res[[cc]][[kk]]$mu[2, , ])
            muList[[cc]][[kk]][[2]] <- t(res[[cc]][[kk]]$mu[1, , ])
            rhoList[[cc]][[kk]][[1]] <- numeric(1000)
            rhoList[[cc]][[kk]][[2]] <- numeric(1000)
            for (ss in 1:1000) {
                rhoList[[cc]][[kk]][[1]][ss] <- cov2cor(res[[cc]][[kk]]$cov[2, , , ss])[1, 2]
                rhoList[[cc]][[kk]][[2]][ss] <- cov2cor(res[[cc]][[kk]]$cov[1, , , ss])[1, 2]
            }
        }
    }
}

saveRDS(rhoList, "~/wasp/mixtures/result/subRhoList10.rds")
saveRDS(muList, "~/wasp/mixtures/result/subMuList10.rds")

for (cc in 1:10) {
    for (kk in 1:nsub) {
        write.table(muList[[cc]][[kk]][[1]], file = paste0("csv/samp_cv_", cc, "_nsub_", kk, "_k10_mu1.csv"),
                    sep = ",", row.names = FALSE, col.names = FALSE)
        write.table(muList[[cc]][[kk]][[2]], file = paste0("csv/samp_cv_", cc, "_nsub_", kk, "_k10_mu2.csv"),
                    sep = ",", row.names = FALSE, col.names = FALSE)
    }
}

library(mvtnorm)
xx <- seq(0, 10, length = 500)
densSampK10 <- rep(list(rep(list(matrix(NA, 1000, 500)), 10)), 10)
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    for (kk in 1:nsub) {
        for (ii in 1:1000) {
            cat("ii: ", ii, "\n")
            for (gg in seq_along(xx)) {
                yy <- c(xx[gg], xx[gg])
                dens1 <- dmvnorm(yy, mean = res[[cc]][[kk]]$mu[1, , ii], sigma = res[[cc]][[kk]]$cov[1, , , ii])
                dens2 <- dmvnorm(yy, mean = res[[cc]][[kk]]$mu[2, , ii], sigma = res[[cc]][[kk]]$cov[2, , , ii])
                densSampK10[[cc]][[kk]][ii, gg] <- res[[cc]][[kk]]$prob[ii, 1] * dens1 + res[[cc]][[kk]]$prob[ii, 2] * dens2
            }
        }
    }
}

saveRDS(densSampK10, "/Shared/ssrivastva/wasp/mixtures/result/waspDensList_k10.rds")

## vb

library(mvtnorm)
xx <- seq(0, 10, length = 500)
vbDens <- list()
vbRho1 <- list()
vbRho2 <- list()
for (cc in 1:10) {
    cat("cc: ", cc, "\n")
    dat <- readRDS(paste0("vb/vb_cv_", cc, ".rds"))
    vbRho1[[cc]] <- numeric(1000)
    vbRho2[[cc]] <- numeric(1000)
    vbDens[[cc]] <- matrix(NA, 1000, 500)
    for (ii in 1:1000) {
        cat("ii: ", ii, "\n")
        infMat1 <- rwish(v = dat$v[1], dat$W[[1]])
        sig1 <- chol2inv(chol(infMat1))
        mu1 <- as.numeric(crossprod(chol(sig1) / sqrt(dat$beta[1]), rnorm(2)) + dat$mu[[1]])
        infMat2 <- rwish(v = dat$v[2], dat$W[[2]])
        sig2 <- chol2inv(chol(infMat2))
        mu2 <- as.numeric(crossprod(chol(sig2) / sqrt(dat$beta[2]), rnorm(2)) + dat$mu[[2]])
        vbRho1[[cc]][ii] <- cov2cor(sig2)[1, 2] # labels are flipped in VB results (2 is 1)
        vbRho2[[cc]][ii] <- cov2cor(sig1)[1, 2] # labels are flipped in VB results (1 is 2)
        pp <- as.numeric(rdirichlet(1, dat$alpha))
        for (gg in seq_along(xx)) {
            yy <- c(xx[gg], xx[gg])
            dens1 <- dmvnorm(yy, mean = mu1, sigma = sig1)
            dens2 <- dmvnorm(yy, mean = mu2, sigma = sig2)
            vbDens[[cc]][ii, gg] <- pp[1] * dens1 + pp[2] * dens2
        }
    }
}

saveRDS(vbRho1, "/Shared/ssrivastva/wasp/mixtures/result/vbRho1.rds")
saveRDS(vbRho2, "/Shared/ssrivastva/wasp/mixtures/result/vbRho2.rds")
saveRDS(vbDens, "/Shared/ssrivastva/wasp/mixtures/result/vbDensList.rds")
