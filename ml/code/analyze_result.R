rm(list=ls())

setwd("/Shared/ssrivastva/wasp/ml/result/")

library(KernSmooth)
library(matrixStats)
library(MCMCpack)
library(RColorBrewer)
library(matrixStats)

colors <- brewer.pal(6, "Set1")

mcmcFix <- list()
scotFix <- list()
xingFix <- list()
waspFix <- list()
vbFix <- list()
bbvbFix <- list()

mcmcRan <- list()
scotRan <- list()
xingRan <- list()
waspRan <- list()
vbRan <- list()
bbvbRan <- list()

mcmcTime <- list()
xingTime <- list()
scotTime <- list()
waspTime <- list()
vbTime <- list()
bbvbTime <- list()

meanMap <- 1:6
covMap <- c(7:12, 14:18, 21:24, 28:30, 35:36, 42)
for (cc in 1:10) {
    mcmcFix[[cc]] <- readRDS(paste0("full/full_res_", cc, ".rds"))$samples[ , meanMap]
    mcmcRan[[cc]] <- readRDS(paste0("full/full_res_", cc, ".rds"))$samples[ , covMap]
    mcmcTime[[cc]] <- readRDS(paste0("full/full_res_", cc, ".rds"))$time[3]
}

meanList <- list()
covList <- list()
for (cc in 1:10) {
    waspTime[[cc]] <- list()
    meanList[[cc]] <- list()
    covList[[cc]] <- list()
    for (kk in 1:10) {
        meanList[[cc]][[kk]] <- readRDS(paste0("wasp/samp/wasp_cv_", cc, "_sub_", kk, "_k10.rds"))$samples[ , meanMap]
        covList[[cc]][[kk]] <- readRDS(paste0("wasp/samp/wasp_cv_", cc, "_sub_", kk, "_k10.rds"))$samples[ , covMap]
        waspTime[[cc]][[kk]] <- readRDS(paste0("wasp/samp/wasp_cv_", cc, "_sub_", kk, "_k10.rds"))$time[3]
    }
}

waspRan <- list()
waspFix <- list()
for (cc in 1:10) {
    waspRan[[cc]] <- list()
    waspFix[[cc]] <- list()
    for (dd in 1:6) {
        waspFix[[cc]][[dd]] <- rowMeans(do.call(cbind, lapply(meanList[[cc]], function(x) quantile(x[ , dd], probs = seq(0, 1, by = 0.0001)))))
    }
    for (dd in 1:21) {
        waspRan[[cc]][[dd]] <- rowMeans(do.call(cbind, lapply(covList[[cc]], function(x) quantile(x[ , dd], probs = seq(0, 1, by = 0.0001)))))
    }
}

for (cc in 1:10) {
    dat <- readRDS(paste0("xing/marg/xing_fix_ran_cv_", cc, "_k10.rds"))
    xingFix[[cc]] <- dat$fix
    xingRan[[cc]] <- dat$ran
    xingTime[[cc]] <- dat$time[3]
}

for (cc in 1:10) {
    dat <- readRDS(paste0("cons/marg/cons_fix_ran_cv_", cc, "_k10.rds"))
    scotFix[[cc]] <- dat$fix
    scotRan[[cc]] <- dat$ran
    scotTime[[cc]] <- dat$time[3]
}

for (cc in 1:10) {
    dat <- readRDS(paste0("vb/vb_res_", cc, ".rds"))
    cmat <- matrix(NA, 1000, sum(lower.tri(dat$cov$scale, diag = TRUE)))
    for (ll in 1:1000) {
        tmp <- riwish(dat$cov$df, dat$cov$scale)
        cmat[ll, ] <- tmp[lower.tri(tmp, diag = TRUE)]
    }
    vbRan[[cc]] <- cmat
    vbFix[[cc]] <- t(crossprod(chol(dat$coefs$cov), matrix(rnorm(1000 * ncol(dat$coefs$cov)), ncol(dat$coefs$cov), 1000)) + dat$coefs$mu)
    vbTime[[cc]] <- dat$time[3]
}

for (cc in 1:10) {
    dat <- readRDS(paste0("bbvb/bbvb_res_", cc, ".rds"))
    lst <- dat[[1]]@sim$samples[[1]]
    bs <- grep("fixef|covRanef", names(lst))
    sampdf <- do.call(cbind, lst[bs])
    bbvbFix[[cc]] <- sampdf[ , meanMap]
    bbvbRan[[cc]] <- sampdf[ , covMap]
    vtmp <- bbvbRan[[cc]][ , corrMap]
    idx1 <- rep(1:6, times = rev(1:6) - 1)
    idx2 <- c(unlist(lapply(2:5, function(x) seq(x, 6))), 6)
    svmat <- sqrt(vtmp[ , idx1] * vtmp[ , idx2])
    bbvbCorr[[cc]] <- bbvbRan[[cc]][ , -corrMap] / svmat
    bbvbTime[[cc]] <- dat$time[3]
}

resRan <- vector("list", 6)
names(resRan) <- c("MCMC", "Xing", "Scot", "WASP", "VB", "bbvb")
for (cc in 1:10) {
    resRan[["MCMC"]][[cc]] <- list()
    resRan[["Xing"]][[cc]] <- list()
    resRan[["Scot"]][[cc]] <- list()
    resRan[["WASP"]][[cc]] <- list()
    resRan[["VB"]][[cc]] <- list()
    resRan[["bbvb"]][[cc]] <- list()
    for (dd in 1:ncol(mcmcRan[[cc]])) {
        rr <- range(c(mcmcRan[[cc]][ , dd],
                      xingRan[[cc]][ , dd],
                      scotRan[[cc]][ , dd],
                      vbRan[[cc]][ , dd],
                      waspRan[[cc]][[dd]],
                      bbvbRan[[cc]][ , dd]))
        bw1 <- dpik(mcmcRan[[cc]][ , dd], range.x = rr, gridsize = 1000)
        bw2 <- dpik(xingRan[[cc]][ , dd], range.x = rr, gridsize = 1000)
        bw3 <- dpik(scotRan[[cc]][ , dd], range.x = rr, gridsize = 1000)
        bw4 <- dpik(waspRan[[cc]][[dd]], range.x = rr, gridsize = 1000)
        bw5 <- dpik(vbRan[[cc]][ , dd], range.x = rr, gridsize = 1000)
        bw6 <- dpik(bbvbRan[[cc]][ , dd], range.x = rr, gridsize = 1000)
        dens1 <- bkde(mcmcRan[[cc]][ , dd], bandwidth = bw1, range.x = rr, gridsize = 1000)
        dens2 <- bkde(xingRan[[cc]][ , dd], bandwidth = bw2, range.x = rr, gridsize = 1000)
        dens3 <- bkde(scotRan[[cc]][ , dd], bandwidth = bw3, range.x = rr, gridsize = 1000)
        dens4 <- bkde(waspRan[[cc]][[dd]], bandwidth = bw4, range.x = rr, gridsize = 1000)
        dens5 <- bkde(vbRan[[cc]][ , dd], bandwidth = bw5, range.x = rr, gridsize = 1000)
        dens6 <- bkde(bbvbRan[[cc]][ , dd], bandwidth = bw6, range.x = rr, gridsize = 1000)
        resRan[["MCMC"]][[cc]][[dd]] <- dens1
        resRan[["Xing"]][[cc]][[dd]] <- dens2
        resRan[["Scot"]][[cc]][[dd]] <- dens3
        resRan[["WASP"]][[cc]][[dd]] <- dens4
        resRan[["VB"]][[cc]][[dd]] <- dens5
        resRan[["bbvb"]][[cc]][[dd]] <- dens6
    }
}

resFix <- vector("list", 6)
names(resFix) <- c("MCMC", "Xing", "Scot", "WASP", "VB", "bbvb")
for (cc in 1:10) {
    resFix[["MCMC"]][[cc]] <- list()
    resFix[["Xing"]][[cc]] <- list()
    resFix[["Scot"]][[cc]] <- list()
    resFix[["WASP"]][[cc]] <- list()
    resFix[["VB"]][[cc]] <- list()
    resFix[["bbvb"]][[cc]] <- list()
    for (dd in 1:ncol(mcmcFix[[cc]])) {
        rr <- range(c(mcmcFix[[cc]][ , dd],
                      xingFix[[cc]][ , dd],
                      scotFix[[cc]][ , dd],
                      vbFix[[cc]][ , dd],
                      waspFix[[cc]][[dd]],
                      bbvbFix[[cc]][ , dd]))
        bw1 <- dpik(mcmcFix[[cc]][ , dd], range.x = rr)
        bw2 <- dpik(xingFix[[cc]][ , dd], range.x = rr)
        bw3 <- dpik(scotFix[[cc]][ , dd], range.x = rr)
        bw4 <- dpik(waspFix[[cc]][[dd]], range.x = rr)
        bw5 <- dpik(vbFix[[cc]][ , dd], range.x = rr)
        bw6 <- dpik(bbvbFix[[cc]][ , dd], range.x = rr)
        dens1 <- bkde(mcmcFix[[cc]][ , dd], bandwidth = bw1, range.x = rr)
        dens2 <- bkde(xingFix[[cc]][ , dd], bandwidth = bw2, range.x = rr)
        dens3 <- bkde(scotFix[[cc]][ , dd], bandwidth = bw3, range.x = rr)
        dens4 <- bkde(waspFix[[cc]][[dd]], bandwidth = bw4, range.x = rr)
        dens5 <- bkde(vbFix[[cc]][ , dd], bandwidth = bw5, range.x = rr)
        dens6 <- bkde(bbvbFix[[cc]][ , dd], bandwidth = bw6, range.x = rr)
        resFix[["MCMC"]][[cc]][[dd]] <- dens1
        resFix[["Xing"]][[cc]][[dd]] <- dens2
        resFix[["Scot"]][[cc]][[dd]] <- dens3
        resFix[["WASP"]][[cc]][[dd]] <- dens4
        resFix[["VB"]][[cc]][[dd]] <- dens5
        resFix[["bbvb"]][[cc]][[dd]] <- dens6
    }
}

accRan <- array(NA,
                dim = c(10, ncol = ncol(mcmcRan[[1]]), 5),
                dimnames = list(paste0("cv", 1:10),
                                paste0("dim", 1:ncol(mcmcRan[[1]])),
                                c("xing", "scot", "wasp", "vb", "bbvb")
                                )
                )
for (cc in 1:10) {
    for (dd in 1:ncol(mcmcRan[[cc]])) {
        accRan[cc, dd, 1] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[dd]]$y  - resRan[["Xing"]][[cc]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accRan[cc, dd, 2] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[dd]]$y  - resRan[["Scot"]][[cc]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accRan[cc, dd, 3] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[dd]]$y  - resRan[["WASP"]][[cc]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accRan[cc, dd, 4] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[dd]]$y  - resRan[["VB"]][[cc]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accRan[cc, dd, 5] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[dd]]$y  - resRan[["bbvb"]][[cc]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
    }
}

accFix <- array(NA,
                dim = c(10, ncol = ncol(mcmcFix[[1]]), 5),
                dimnames = list(paste0("cv", 1:10),
                                paste0("dim", 1:ncol(mcmcFix[[1]])),
                                c("xing", "scot", "wasp", "vb", "bbvb")
                                )
                )
for (cc in 1:10) {
    for (dd in 1:ncol(mcmcFix[[cc]])) {
        accFix[cc, dd, 1] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[dd]]$y  - resFix[["Xing"]][[cc]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accFix[cc, dd, 2] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[dd]]$y  - resFix[["Scot"]][[cc]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accFix[cc, dd, 3] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[dd]]$y  - resFix[["WASP"]][[cc]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accFix[cc, dd, 4] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[dd]]$y  - resFix[["VB"]][[cc]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
        accFix[cc, dd, 5] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[dd]]$y  - resFix[["bbvb"]][[cc]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[dd]]$x)[1]) / 2)
    }
}

fixTblMn <- format(round(rbind(colMeans(accFix[ , , "scot"]), colMeans(accFix[ , , "xing"]),
                               colMeans(accFix[ , , "vb"]), colMeans(accFix[ , , "bbvb"]), colMeans(accFix[ , , "wasp"])), 2), nsmall = 2)
fixTblSd <- format(round(rbind(colSds(accFix[ , , "scot"]), colSds(accFix[ , , "xing"]),
                               colSds(accFix[ , , "vb"]), colSds(accFix[ , , "bbvb"]), colSds(accFix[ , , "wasp"])), 2), nsmall = 2)
fixTbl <- matrix(paste0(fixTblMn, " (", fixTblSd, ")"), 5, ncol(accFix[ , , "scot"]))

ranTblMn <- format(round(rbind(colMeans(accRan[ , , "scot"]), colMeans(accRan[ , , "xing"]),
                               colMeans(accRan[ , , "vb"]), colMeans(accRan[ , , "bbvb"]),
                               colMeans(accRan[ , , "wasp"])), 2), nsmall = 2)
ranTblSd <- format(round(rbind(colSds(accRan[ , , "scot"]), colSds(accRan[ , , "xing"]),
                               colSds(accRan[ , , "vb"]), colSds(accRan[ , , "bbvb"]),
                               colSds(accRan[ , , "wasp"])), 2), nsmall = 2)
ranTbl <- matrix(paste0(ranTblMn, " (", ranTblSd, ")"), 5, ncol(accRan[ , , "scot"]))

rownames(fixTbl) <- rownames(ranTbl) <- c("CMC", "SDP", "VB", "BBVB", "WASP")

xtable::xtable(fixTbl)
xtable::xtable(ranTbl[ , c(1, 7, 12, 16, 19, 21)])
cidx <- (1:21)[-c(1, 7, 12, 16, 19, 21)]
xtable::xtable(ranTbl[ , cidx[1:8]])
xtable::xtable(ranTbl[ , cidx[-(1:8)]])

xtable::xtable(corrTbl[ , 1:8])
xtable::xtable(corrTbl[ , -(1:8)])

#### 2d plots

mcmcCovMat <- list()
vbCovMat <- list()
scotCovMat <- list()
xingCovMat <- list()
waspCovMat <- list()
bbvbCovMat <- list()

cov2d <- cbind(rep(2, 4), 3:6)
for (cc in 1:10) {
    mcmcCovMat[[cc]] <- list()
    vbCovMat[[cc]] <- list()
    bbvbCovMat[[cc]] <- list()
    for (dd in 1:4) {
        mcmcCovMat[[cc]][[dd]] <- mcmcRan[[cc]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
        bbvbCovMat[[cc]][[dd]] <- bbvbRan[[cc]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
        vbCovMat[[cc]][[dd]] <- vbRan[[cc]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
    }
}

for (cc in 1:10) {
    dat <- readRDS(paste0("xing/joint/xing_cov_cv_", cc, "_k10.rds"))
    xingCovMat[[cc]] <- dat$cov
    xingTime[[cc]] <- dat$time
}

for (cc in 1:10) {
    dat <- readRDS(paste0("cons/joint/cons_cov_cv_", cc, "_k10.rds"))
    scotCovMat[[cc]] <- dat$cov
    scotTime[[cc]] <- dat$time
}

for (cc in 1:10) {
    waspCovMat[[cc]] <- list()
    for (dd in 1:4) {
        dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_d_", dd, "_k10.csv"), header = FALSE)
        pr <- as.numeric(dat[ , 3])
        pr[pr < 1e-10] <- 0.0
        waspCovMat[[cc]][[dd]] <- dat[sample(1:nrow(dat), 10000, replace = TRUE, prob = pr), 1:2]
    }
}

resCovMat <- vector("list", 6)
names(resCovMat) <- c("MCMC", "Xing", "Scot", "WASP", "VB", "bbvb")
for (cc in 1:10) {
    resCovMat[["MCMC"]][[cc]] <- list()
    resCovMat[["Xing"]][[cc]] <- list()
    resCovMat[["Scot"]][[cc]] <- list()
    resCovMat[["WASP"]][[cc]] <- list()
    resCovMat[["VB"]][[cc]] <- list()
    resCovMat[["bbvb"]][[cc]] <- list()
    for (dd in 1:4) {
        dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_d_", dd, "_k10.csv"), header = FALSE)
        bw1 <- max(diff(as.numeric(dat[ , 1]))) * 2
        bw2 <- max(diff(as.numeric(dat[ , 2]))) * 2

        rr1 <- range(c(mcmcCovMat[[cc]][[dd]][ , 1],
                       xingCovMat[[cc]][[dd]][ , 1],
                       scotCovMat[[cc]][[dd]][ , 1],
                       waspCovMat[[cc]][[dd]][ , 1],
                       vbCovMat[[cc]][[dd]][ , 1],
                       bbvbCovMat[[cc]][[dd]][ , 1]
                       ))
        rr2 <- range(c(mcmcCovMat[[cc]][[dd]][ , 2],
                       xingCovMat[[cc]][[dd]][ , 2],
                       scotCovMat[[cc]][[dd]][ , 2],
                       waspCovMat[[cc]][[dd]][ , 2],
                       vbCovMat[[cc]][[dd]][ , 2],
                       bbvbCovMat[[cc]][[dd]][ , 2]
                       ))

        bw12 <- bw22 <- bw32 <- bw42 <- bw52 <- bw62 <- bw2
        bw11 <- bw21 <- bw31 <- bw41 <- bw51 <- bw61 <- bw1
        dens1 <- bkde2D(mcmcCovMat[[cc]][[dd]], bandwidth = c(bw11, bw12), range.x = list(rr1, rr2))
        dens2 <- bkde2D(xingCovMat[[cc]][[dd]], bandwidth = c(bw21, bw22), range.x = list(rr1, rr2))
        dens3 <- bkde2D(scotCovMat[[cc]][[dd]], bandwidth = c(bw31, bw32), range.x = list(rr1, rr2))
        dens4 <- bkde2D(waspCovMat[[cc]][[dd]], bandwidth = c(bw41, bw42), range.x = list(rr1, rr2))
        dens5 <- bkde2D(vbCovMat[[cc]][[dd]], bandwidth = c(bw51, bw52), range.x = list(rr1, rr2))
        dens6 <- bkde2D(bbvbCovMat[[cc]][[dd]], bandwidth = c(bw61, bw62), range.x = list(rr1, rr2))
        resCovMat[["MCMC"]][[cc]][[dd]] <- dens1
        resCovMat[["Xing"]][[cc]][[dd]] <- dens2
        resCovMat[["Scot"]][[cc]][[dd]] <- dens3
        resCovMat[["WASP"]][[cc]][[dd]] <- dens4
        resCovMat[["VB"]][[cc]][[dd]] <- dens5
        resCovMat[["bbvb"]][[cc]][[dd]] <- dens6
    }
}

accCovMat <- array(NA,
                   dim = c(10, 4, 5),
                   dimnames = list(paste0("cv", 1:10),
                                   paste(colnames(mcmcRan[[1]])[cov2d[ , 1]], colnames(mcmcRan[[1]])[cov2d[ , 2]], sep = ","),
                                   c("xing", "scot", "wasp", "vb", "bbvb")
                                   ))
for (cc in 1:10) {
    for (dd in 1:4) {
        accCovMat[cc, dd, 1] <- 1 - sum(abs(resCovMat[["MCMC"]][[cc]][[dd]]$fhat - resCovMat[["Xing"]][[cc]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x2)[1]) / 2
        accCovMat[cc, dd, 2] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[dd]]$fhat - resCovMat[["Scot"]][[cc]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x2)[1]) / 2
        accCovMat[cc, dd, 3] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[dd]]$fhat - resCovMat[["WASP"]][[cc]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x2)[1]) / 2
        accCovMat[cc, dd, 4] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[dd]]$fhat - resCovMat[["VB"]][[cc]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x2)[1]) / 2
        accCovMat[cc, dd, 5] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[dd]]$fhat - resCovMat[["bbvb"]][[cc]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[dd]]$x2)[1]) / 2
    }
}

cov2dTblMn <- format(round(rbind(colMeans(accCovMat[ , , "scot"]), colMeans(accCovMat[ , , "xing"]),
                                  colMeans(accCovMat[ , , "vb"]), colMeans(accCovMat[ , , "bbvb"]),  colMeans(accCovMat[ , , "wasp"])), 2), nsmall = 2)
cov2dTblSd <- format(round(rbind(colSds(accCovMat[ , , "scot"]), colSds(accCovMat[ , , "xing"]),
                                  colSds(accCovMat[ , , "vb"]), colSds(accCovMat[ , , "bbvb"]),         colSds(accCovMat[ , , "wasp"])), 2), nsmall = 2)

cov2dTbl <- matrix(paste0(cov2dTblMn, " (", cov2dTblSd, ")"), 5, ncol(accCovMat[ , , "scot"]))
rownames(cov2dTbl) <- c("CMC", "SDP", "VB", "BBVB", "WASP")
xtable::xtable(cov2dTbl)

rnames <- list(c("12", "13"), c("12", "14"), c("12", "15"), c("12", "16"))

cc <- 1
pdf("~/wasp/ml/result/img/cov_2d.pdf", width = 40, height = 10)
par(mfrow = c(1, 4))
par(cex = 1)
par(mar = c(2, 4, 0, 2), oma = c(1, 1, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (dd in 1:4) {
    contour(resCovMat[["MCMC"]][[cc]][[dd]]$x1, resCovMat[["MCMC"]][[cc]][[dd]]$x2,
            resCovMat[["MCMC"]][[cc]][[dd]]$fhat,
            xlim = range(c(resCovMat[["MCMC"]][[cc]][[dd]]$x1, resCovMat[["CMC"]][[cc]][[dd]]$x1,
                           resCovMat[["SDP"]][[cc]][[dd]]$x1, resCovMat[["VB"]][[cc]][[dd]]$x1,
                           resCovMat[["WASP"]][[cc]][[dd]]$x1)),
            ylim = range(c(resCovMat[["MCMC"]][[cc]][[dd]]$x2, resCovMat[["CMC"]][[cc]][[dd]]$x2,
                           resCovMat[["SDP"]][[cc]][[dd]]$x2, resCovMat[["VB"]][[cc]][[dd]]$x2,
                           resCovMat[["WASP"]][[cc]][[dd]]$x2)) + c(-0.01, 0.01),
            lwd = 6, axes = FALSE, col = colors[1]
            )
    xxlab <- axTicks(1)
    yylab <- axTicks(2)
    contour(resCovMat[["Scot"]][[cc]][[dd]]$x1, resCovMat[["Scot"]][[cc]][[dd]]$x2,
            resCovMat[["Scot"]][[cc]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[2])
    contour(resCovMat[["Xing"]][[cc]][[dd]]$x1, resCovMat[["Xing"]][[cc]][[dd]]$x2,
            resCovMat[["Xing"]][[cc]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[3])
    contour(resCovMat[["VB"]][[cc]][[dd]]$x1, resCovMat[["VB"]][[cc]][[dd]]$x2,
            resCovMat[["VB"]][[cc]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[4])
    contour(resCovMat[["WASP"]][[cc]][[dd]]$x1, resCovMat[["WASP"]][[cc]][[dd]]$x2,
            resCovMat[["WASP"]][[cc]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[5])
    contour(resCovMat[["bbvb"]][[cc]][[dd]]$x1, resCovMat[["bbvb"]][[cc]][[dd]]$x2,
            resCovMat[["bbvb"]][[cc]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[6])
    axis(side = 2, tck = -0.01, lwd = 4, cex = 2, labels = NA)
    mtext(format(round(yylab, 2)), at = yylab, side = 2, line = 0.5, cex = 2.5, las = 2)
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
    grid(lwd = 2)
    box(col = "grey40", lwd = 4)
    mtext(bquote(sigma[.(rnames[[dd]][1])]), side = 3, line = -4, adj = 0.1, cex = 4)
    mtext(", ", side = 3, line = -3.5, adj = 0.22, cex = 4)
    mtext(bquote(sigma[.(rnames[[dd]][2])]), side = 3, line = -4, adj = 0.28, cex = 4)
    if (dd == 3) {
        legend("bottom", c("ADVI", "CMC", "MCMC", "SA", "SDP", "WASP"), lty = c("solid"),
               col = colors[c(6, 2, 1, 4, 3, 5)], lwd = 8, bty = "n", cex = 2.5, ncol = 3)
    }
}
dev.off()

for (cc in 1:10) {
    for (kk in 1:10) {
        dat <- read.csv(paste0("wasp/joint/cov_2d_times_cv_", cc, "_k10.csv"), header = FALSE)
        waspTime[[cc]][[kk]] <- waspTime[[cc]][[kk]] + mean(unlist(dat))
    }
}

timeMat <- cbind.data.frame("mcmc" = unlist(mcmcTime),
                            "bbvb" = unlist(bbvbTime),
                            "vb" = unlist(vbTime),
                            "scot" = unlist(lapply(scotTime, mean)),
                            "xing" = unlist(lapply(xingTime, mean)),
                            "wasp" = unlist(lapply(waspTime, function(x) mean(unlist(x))))
                            )

pdf("~/wasp/ml/result/img/ml_time.pdf", 8, 8)
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(4, 4.5, 0.2, 0.2))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
boxplot(log10(timeMat), ylab = NA, axes = FALSE, lwd = 2, ylim = c(2, 5.1))
box(col = "grey40", lwd = 4)
grid(lwd = 3)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("MCMC", "ADVI", "SA", "CMC", "SDP", "WASP"), at = 1:6, side = 1, line = 1, cex = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(seq(2, 5, by = 0.5)), at = seq(2, 5, by = 0.5), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression(log[10] ~ " Seconds"), side = 2, outer = TRUE, cex = 2, line = 2.7)
dev.off()
