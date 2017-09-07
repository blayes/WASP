rm(list = ls())
setwd("/Shared/ssrivastva/wasp/mixtures/result/")

library(KernSmooth)
library(matrixStats)
library(xtable)

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

wasp10Mu1 <- list()
for (cc in 1:10) {
    dat <- read.csv(paste0("sub10/wasp_cv_", cc, "_mu1_k10.csv"), header = FALSE)
    colnames(dat) <- c("dim1", "dim2", "prob")
    dat$prob[dat$prob < 1e-10] <- 0
    wasp10Mu1[[cc]] <- dat[sample(1:nrow(dat), size = 10000, prob = dat$prob, replace = TRUE), 1:2]
}

wasp10Mu2 <- list()
for (cc in 1:10) {
    dat <- read.csv(paste0("sub10/wasp_cv_", cc, "_mu2_k10.csv"), header = FALSE)
    colnames(dat) <- c("dim1", "dim2", "prob")
    dat$prob[dat$prob < 1e-10] <- 0
    wasp10Mu2[[cc]] <- dat[sample(1:nrow(dat), size = 10000, prob = dat$prob, replace = TRUE), 1:2]
}

wasp5Mu1 <- list()
for (cc in 1:10) {
    dat <- read.csv(paste0("sub5/wasp_cv_", cc, "_mu1_k5.csv"), header = FALSE)
    colnames(dat) <- c("dim1", "dim2", "prob")
    dat$prob[dat$prob < 1e-10] <- 0
    wasp5Mu1[[cc]] <- dat[sample(1:nrow(dat), size = 10000, prob = dat$prob, replace = TRUE), 1:2]
}

wasp5Mu2 <- list()
for (cc in 1:10) {
    dat <- read.csv(paste0("sub5/wasp_cv_", cc, "_mu2_k5.csv"), header = FALSE)
    colnames(dat) <- c("dim1", "dim2", "prob")
    dat$prob[dat$prob < 1e-10] <- 0
    wasp5Mu2[[cc]] <- dat[sample(1:nrow(dat), size = 10000, prob = dat$prob, replace = TRUE), 1:2]
}

accMu <- array(0.0, dim = c(10, 2, 2))
dimnames(accMu) <- list(c(paste0("cv", 1:10)),
                         c("k5", "k10"),
                         c("mu1", "mu2")
                        )

for (cc in 1:10) {
        waspDens51 <- bkde2D(wasp5Mu1[[cc]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                            range.x =  list(c(0.95, 1.05), c(1.95, 2.05)))
        waspDens52 <- bkde2D(wasp5Mu2[[cc]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                             range.x =  list(c(6.95, 7.05), c(7.95, 8.05)))
        waspDens101 <- bkde2D(wasp10Mu1[[cc]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                            range.x =  list(c(0.95, 1.05), c(1.95, 2.05)))
        waspDens102 <- bkde2D(wasp10Mu2[[cc]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                            range.x =  list(c(6.95, 7.05), c(7.95, 8.05)))
        fullDens1 <- bkde2D(muList[[cc]][[1]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                            range.x =  list(c(0.95, 1.05), c(1.95, 2.05)))
        fullDens2 <- bkde2D(muList[[cc]][[2]], bandwidth = c(0.01, 0.01), gridsize = c(500, 500),
                            range.x =  list(c(6.95, 7.05), c(7.95, 8.05)))
        accMu[cc, 1, 1] <- 1- sum(abs(fullDens1$fhat - waspDens51$fhat) * diff(fullDens1$x1)[1] * diff(fullDens1$x2)[1]) / 2
        accMu[cc, 1, 2] <- 1- sum(abs(fullDens2$fhat - waspDens52$fhat) * diff(fullDens2$x1)[1] * diff(fullDens2$x2)[1]) / 2
        accMu[cc, 2, 1] <- 1- sum(abs(fullDens1$fhat - waspDens101$fhat) * diff(fullDens1$x1)[1] * diff(fullDens1$x2)[1]) / 2
        accMu[cc, 2, 2] <- 1- sum(abs(fullDens2$fhat - waspDens102$fhat) * diff(fullDens2$x1)[1] * diff(fullDens2$x2)[1]) / 2
}

saveRDS(accMu, "~/wasp/mixtures/result/accMu.rds")

## see create_samples.R file for their description and structure.
fullRho <- readRDS("~/wasp/mixtures/result/corrList.rds")
rho5 <- readRDS("~/wasp/mixtures/result/subRhoList5.rds")
rho10 <- readRDS( "~/wasp/mixtures/result/subRhoList10.rds")
vbRho1 <- readRDS( "/Shared/ssrivastva/wasp/mixtures/result/vbRho1.rds")
vbRho2 <- readRDS( "/Shared/ssrivastva/wasp/mixtures/result/vbRho2.rds")

scotRho <- list()
for (cc in 1:10) {
    scotRho[[cc]] <- list()
    dat1 <- readRDS(paste0("cons/cons_rho_cv_", cc, "_k5.rds"))
    dat2 <- readRDS(paste0("cons/cons_rho_cv_", cc, "_k10.rds"))
    scotRho[[cc]] <- list("k5" = dat1$dens, "k10" = dat2$dens)
}

xingRho <- list()
for (cc in 1:10) {
    dat1 <- readRDS(paste0("xing/xing_rho_cv_", cc, "_k5.rds"))
    dat2 <- readRDS(paste0("xing/xing_rho_cv_", cc, "_k10.rds"))
    xingRho[[cc]] <- list("k5" = dat1$dens, "k10" = dat2$dens)
}

waspRho5 <- list(lapply(rho5,
                        function (x) {
                            rowMeans(do.call(cbind, lapply(x, function(y) quantile(y[[1]], probs = seq(0, 1, by = 0.0001)))))
                        }),
                 lapply(rho5,
                        function (x) {
                            rowMeans(do.call(cbind, lapply(x, function(y) quantile(y[[2]], probs = seq(0, 1, by = 0.0001)))))
                        }))

waspRho10 <- list(lapply(rho10,
                        function (x) {
                            rowMeans(do.call(cbind, lapply(x, function(y) quantile(y[[1]], probs = seq(0, 1, by = 0.0001)))))
                        }),
                 lapply(rho10,
                        function (x) {
                            rowMeans(do.call(cbind, lapply(x, function(y) quantile(y[[2]], probs = seq(0, 1, by = 0.0001)))))
                        }))

accRho <- array(NA, dim = c(10, 4, 2, 2))
dimnames(accRho) <- list(c(paste0("cv", 1:10)),
                         c("cons", "xing", "vb", "wasp"),
                         c("rho1", "rho2"),
                         c("k5", "k10")
                        )

for (cc in 1:10) {
    for (rrr in 1:2) {
        rr <- range(c(fullRho[[cc]][[rrr]], waspRho5[[rrr]][[cc]], waspRho10[[rrr]][[cc]],
                      as.numeric(xingRho[[cc]][[1]][ , rrr]), as.numeric(xingRho[[cc]][[2]][ , rrr]),
                      as.numeric(scotRho[[cc]][[1]][ , rrr]), as.numeric(scotRho[[cc]][[2]][ , rrr]),
                      as.numeric(vbRho1[[cc]]), as.numeric(vbRho2[[cc]])
                      ))
        fdens <- bkde(fullRho[[cc]][[rrr]], bandwidth = dpik(fullRho[[cc]][[rrr]]), range.x = rr)
        wdens1 <- bkde(waspRho5[[rrr]][[cc]], bandwidth = dpik(waspRho5[[rrr]][[cc]]), range.x = rr)
        wdens2 <- bkde(waspRho10[[rrr]][[cc]], bandwidth = dpik(waspRho10[[rrr]][[cc]]), range.x = rr)
        xdens1 <- bkde(xingRho[[cc]][[1]][ , rrr], bandwidth = dpik(xingRho[[cc]][[1]][ , rrr]), range.x = rr)
        xdens2 <- bkde(xingRho[[cc]][[2]][ , rrr], bandwidth = dpik(xingRho[[cc]][[2]][ , rrr]), range.x = rr)
        cdens1 <- bkde(scotRho[[cc]][[1]][ , rrr], bandwidth = dpik(scotRho[[cc]][[1]][ , rrr]), range.x = rr)
        cdens2 <- bkde(scotRho[[cc]][[2]][ , rrr], bandwidth = dpik(scotRho[[cc]][[2]][ , rrr]), range.x = rr)
        if (rrr == 1) {
            vdens <- bkde(vbRho1[[cc]], bandwidth = dpik(vbRho1[[cc]]), range.x = rr)
        } else {
            vdens <- bkde(vbRho2[[cc]], bandwidth = dpik(vbRho2[[cc]]), range.x = rr)
        }
        accRho[cc, "wasp", rrr, 1] <- 1- sum(abs(fdens$y - wdens1$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "wasp", rrr, 2] <- 1- sum(abs(fdens$y - wdens2$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "xing", rrr, 1] <- 1- sum(abs(fdens$y - xdens1$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "xing", rrr, 2] <- 1- sum(abs(fdens$y - xdens2$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "cons", rrr, 1] <- 1- sum(abs(fdens$y - cdens1$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "cons", rrr, 2] <- 1- sum(abs(fdens$y - cdens2$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "vb", rrr, 1] <- 1- sum(abs(fdens$y - vdens$y) * diff(fdens$x)[1]) / 2
        accRho[cc, "vb", rrr, 2] <- 1- sum(abs(fdens$y - vdens$y) * diff(fdens$x)[1]) / 2
    }
}

saveRDS(accRho, "~/wasp/mixtures/result/accRho.rds")

rtbl <- rbind(
    c(paste0(format(round(colMeans(accRho[ , "cons", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "cons", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accRho[ , "cons", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "cons", , 2]), 2), nsmall = 2), ")")
      ),
    c(paste0(format(round(colMeans(accRho[ , "xing", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "xing", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accRho[ , "xing", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "xing", , 2]), 2), nsmall = 2), ")")
      ),
    c(paste0(format(round(colMeans(accRho[ , "vb", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "vb", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accRho[ , "vb", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "vb", , 2]), 2), nsmall = 2), ")")
      ),
    c(paste0(format(round(colMeans(accRho[ , "wasp", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "wasp", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accRho[ , "wasp", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accRho[ , "wasp", , 2]), 2), nsmall = 2), ")")
      )
)

rownames(rtbl) <- c("cons", "xing", "vb", "wasp")
colnames(rtbl) <- c("rho1-k5", "rho2-k5", "rho1-k10", "rho2-k10")

xtable(rtbl[, c(1, 3, 2, 4)])
xtable(resTbl)

### f-estimation ###

rm(list = ls())

setwd("/Shared/ssrivastva/wasp/mixtures/result/")

mcmcDens <- readRDS("mcmcDensList.rds")
vbDens <- readRDS("vbDensList.rds")
waspDens10 <- readRDS("waspDensList_k10.rds")
waspDens5 <- readRDS("waspDensList_k5.rds")

waspDens <- list()
for (cc in 1:10) {
    waspDens[[cc]] <- list()
    densMat <- matrix(NA, 2000, 500)
    for (ss in 1:500) {
        densMat[, ss] <- rowMeans(do.call(cbind, lapply(lapply(waspDens5[[cc]], function(x) x[ , ss]),
                                                         function(y) quantile(y, seq(0, 1, length = 2000)))))
    }
    waspDens[[cc]][[1]] <- densMat
    densMat <- matrix(NA, 2000, 500)
    for (ss in 1:500) {
        densMat[, ss] <- rowMeans(do.call(cbind, lapply(lapply(waspDens10[[cc]], function(x) x[ , ss]),
                                                         function(y) quantile(y, seq(0, 1, length = 2000)))))
    }
    waspDens[[cc]][[2]] <- densMat
}

scotDens <- list()
for (cc in 1:10) {
    scotDens[[cc]] <- list()
    dat1 <- readRDS(paste0("cons/cons_dens_cv_", cc, "_k5.rds"))
    dat2 <- readRDS(paste0("cons/cons_dens_cv_", cc, "_k10.rds"))
    scotDens[[cc]] <- list("k5" = dat1$dens, "k10" = dat2$dens)
}

xingDens <- list()
for (cc in 1:10) {
    dat1 <- readRDS(paste0("xing/xing_dens_cv_", cc, "_k5.rds"))
    dat2 <- readRDS(paste0("xing/xing_dens_cv_", cc, "_k10.rds"))
    xingDens[[cc]] <- list("k5" = dat1$dens, "k10" = dat2$dens)
}

accMix <- array(NA, dim = c(10, 3, 4, 2))
dimnames(accMix) <- list(paste0("cv", 1:10),
                         c("cons", "vb", "wasp"),
                         c("2.5", "5", "90", "95"),
                         c("k=5", "k=10")
                         )

xx <- seq(0, 10, length = 500)
for (cc in 1:10) {
    mf <- colQuantiles(mcmcDens[[cc]], probs = c(0.025, 0.05, 0.95, 0.975))
    wf5 <- colQuantiles(waspDens[[cc]][[1]], probs = c(0.025, 0.05, 0.95, 0.975))
    wf10 <- colQuantiles(waspDens[[cc]][[2]], probs = c(0.025, 0.05, 0.95, 0.975))
    cf5 <- colQuantiles(scotDens[[cc]][[1]], probs = c(0.025, 0.05, 0.95, 0.975))
    cf10 <- colQuantiles(scotDens[[cc]][[2]], probs = c(0.025, 0.05, 0.95, 0.975))
    vf5 <- colQuantiles(vbDens[[cc]], probs = c(0.025, 0.05, 0.95, 0.975))
    vf10 <- colQuantiles(vbDens[[cc]], probs = c(0.025, 0.05, 0.95, 0.975))
    for (dd in 1:4) {
        accMix[cc, "cons", dd, 1] <- 1 - sum(diff(xx)[1] * abs(mf - cf5[ , dd])) / 2
        accMix[cc, "cons", dd, 2] <- 1 - sum(diff(xx)[1] * abs(mf - cf10[ , dd])) / 2
        accMix[cc, "vb", dd, 1] <- 1 - sum(diff(xx)[1] * abs(mf - vf5[ , dd])) / 2
        accMix[cc, "vb", dd, 2] <- 1 - sum(diff(xx)[1] * abs(mf - vf10[ , dd])) / 2
        accMix[cc, "wasp", dd, 1] <- 1 - sum(diff(xx)[1] * abs(mf - wf5[ , dd])) / 2
        accMix[cc, "wasp", dd, 2] <- 1 - sum(diff(xx)[1] * abs(mf - wf10[ , dd])) / 2
    }
}

rbind(
    c(paste0(format(round(colMeans(accMix[ , "cons", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "cons", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accMix[ , "cons", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "cons", , 2]), 2), nsmall = 2), ")")
      ),
    c(paste0(format(round(colMeans(accMix[ , "vb", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "vb", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accMix[ , "vb", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "vb", , 2]), 2), nsmall = 2), ")")
      ),
    c(paste0(format(round(colMeans(accMix[ , "wasp", , 1]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "wasp", , 1]), 2), nsmall = 2), ")"),
      paste0(format(round(colMeans(accMix[ , "wasp", , 2]), 2), nsmall = 2),
             " (", format(round(colSds(accMix[ , "wasp", , 2]), 2), nsmall = 2), ")")
      )
)
