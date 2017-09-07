rm(list=ls())

setwd("/Shared/ssrivastva/wasp/lme/result/")

library(KernSmooth)
library(matrixStats)
library(MCMCpack)
library(RColorBrewer)
library(rstan)
library(Rcpp)
colors <- brewer.pal(6, "Set1")

mcmcFix <- list()
vbFix <- list()
scotFix <- list()
xingFix <- list()
waspFixSamp <- list()
bbvbFix <- list()

mcmcRan <- list()
vbRan <- list()
scotRan <- list()
xingRan <- list()
waspRanSamp <- list()
bbvbRan <- list()

mcmcTime <- list()
vbTime <- list()
xingTime <- list()
scotTime <- list()
waspTime <- list()
bbvbTime <- list()

ndim <- 1:2
npart <- c(10, 20)

for (cc in 1:10) {
    mcmcRan[[cc]] <- list()
    mcmcFix[[cc]] <- list()
    mcmcTime[[cc]] <- list()
    for (pp in 1:2) {
        dat <- readRDS(paste0("mcmc/mcmc_lme_cv_", cc, "_p_", ndim[pp], ".rds"))
        cnames <- colnames(dat$samples)
        mcmcFix[[cc]][[pp]] <- dat$samples[ , grep("fixef", cnames)]
        if (pp == 2) {
            mcmcRan[[cc]][[pp]] <- dat$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
        } else {
            mcmcRan[[cc]][[pp]] <- dat$samples[ , c(5:7, 9:10, 13)]
        }
        mcmcTime[[cc]][[pp]] <- dat$time[3]
    }
}

for (cc in 1:10) {
    bbvbRan[[cc]] <- list()
    bbvbFix[[cc]] <- list()
    bbvbTime[[cc]] <- list()
    for (pp in 1:2) {
        dat <- readRDS(paste0("bbvb/bbvb_lme_cv_", cc, "_p_", ndim[pp], ".rds"))
        lst <- dat[[1]]@sim$samples[[1]]
        bs <- grep("fixef|covRanef", names(lst))
        sampdf <- do.call(cbind, lst[bs])
        cnames <- colnames(sampdf)
        bbvbFix[[cc]][[pp]] <- sampdf[ , grep("fixef", cnames)]
        if (pp == 2) {
            bbvbRan[[cc]][[pp]] <- sampdf[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
        } else {
            bbvbRan[[cc]][[pp]] <- sampdf[ , c(5:7, 9:10, 13)]
        }
        bbvbTime[[cc]][[pp]] <- dat$time[3]
    }
}

for (cc in 1:10) {
    vbRan[[cc]] <- list()
    vbFix[[cc]] <- list()
    vbTime[[cc]] <- list()
    for (pp in 1:2) {
        dat <- readRDS(paste0("vb/vb_lme_cv_", cc, "_p_", ndim[pp], ".rds"))
        cmat <- matrix(NA, 1000, sum(lower.tri(dat$cov$scale, diag = TRUE)))
        rmat <- matrix(NA, 1000, sum(lower.tri(dat$cov$scale, diag = FALSE)))
        for (ll in 1:1000) {
            tmp <- riwish(dat$cov$df, dat$cov$scale)
            cmat[ll, ] <- tmp[lower.tri(tmp, diag = TRUE)]
            rmat[ll, ] <- cov2cor(tmp)[lower.tri(tmp, diag = FALSE)]
        }
        vbRan[[cc]][[pp]] <- cmat
        vbFix[[cc]][[pp]] <- t(crossprod(chol(dat$coefs$cov), matrix(rnorm(1000 * ncol(dat$coefs$cov)), ncol(dat$coefs$cov), 1000)) + dat$coefs$mu)
        vbTime[[cc]][[pp]] <- dat$time[3]
    }
}

waspTime10 <- list()
for (cc in 1:10) {
    waspRanSamp[[cc]] <- list()
    waspFixSamp[[cc]] <- list()
    waspTime10[[cc]] <- list()
    for (pp in 1:2) {
        waspFixSamp[[cc]][[pp]] <- list()
        waspRanSamp[[cc]][[pp]] <- list()
        waspTime10[[cc]][[pp]] <- list()
        for (kk in 1:10) {
            cat("loaded: ", paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", ndim[pp], "_k_", kk, "_nsub10.rds"), "\n")
            dat <- readRDS(paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", ndim[pp], "_k_", kk, "_nsub10.rds"))
            waspTime10[[cc]][[pp]][[kk]] <- dat$time[3]
            cnames <- colnames(dat$samples)
            waspFixSamp[[cc]][[pp]][[kk]] <- dat$samples[ , grep("fixef", cnames)]
            if (pp == 2) {
                waspRanSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            } else {
                waspRanSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(5:7, 9:10, 13)]
            }
        }
    }
}

waspRan10 <- list()
waspFix10 <- list()
for (cc in 1:10) {
    waspRan10[[cc]] <- list()
    waspFix10[[cc]] <- list()
    for (pp in 1:2) {
        tmp <- list()
        for (kk in 1:10) {
            tmp[[kk]] <- do.call(cbind, lapply(split(waspRanSamp[[cc]][[pp]][[kk]], col(waspRanSamp[[cc]][[pp]][[kk]])),
                                               function(x) quantile(x, probs = seq(0, 1, length = 1000))))
        }
        waspRan10[[cc]][[pp]] <- matrix(NA, nrow = nrow(tmp[[1]]), ncol = ncol(vbRan[[cc]][[pp]]))
        for (dd in 1:ncol(waspRan10[[cc]][[pp]])) {
            waspRan10[[cc]][[pp]][ , dd] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , dd])))
        }
        tmp <- list()
        for (kk in 1:10) {
            tmp[[kk]] <- do.call(cbind, lapply(split(waspFixSamp[[cc]][[pp]][[kk]], col(waspFixSamp[[cc]][[pp]][[kk]])),
                                               function(x) quantile(x, probs = seq(0, 1, length = 1000))))
        }
        waspFix10[[cc]][[pp]] <- matrix(NA, nrow = nrow(tmp[[1]]), ncol = ncol(vbFix[[cc]][[pp]]))
        for (dd in 1:ncol(waspFix10[[cc]][[pp]])) {
            waspFix10[[cc]][[pp]][ , dd] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , dd])))
        }
    }
}

waspTime20 <- list()
for (cc in 1:10) {
    waspRanSamp[[cc]] <- list()
    waspFixSamp[[cc]] <- list()
    waspTime20[[cc]] <- list()
    for (pp in 1:2) {
        waspFixSamp[[cc]][[pp]] <- list()
        waspRanSamp[[cc]][[pp]] <- list()
        waspTime20[[cc]][[pp]] <- list()
        for (kk in 1:20) {
            cat("loaded: ", paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", ndim[pp], "_k_", kk, "_nsub20.rds"), "\n")
            dat <- readRDS(paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", ndim[pp], "_k_", kk, "_nsub20.rds"))
            waspTime20[[cc]][[pp]][[kk]] <- dat$time[3]
            cnames <- colnames(dat$samples)
            waspFixSamp[[cc]][[pp]][[kk]] <- dat$samples[ , grep("fixef", cnames)]
            if (pp == 2) {
                waspRanSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            } else {
                waspRanSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(5:7, 9:10, 13)]
            }
        }
    }
}

waspRan20 <- list()
waspFix20 <- list()
for (cc in 1:10) {
    waspRan20[[cc]] <- list()
    waspFix20[[cc]] <- list()
    for (pp in 1:2) {
        tmp <- list()
        for (kk in 1:20) {
            tmp[[kk]] <- do.call(cbind, lapply(split(waspRanSamp[[cc]][[pp]][[kk]], col(waspRanSamp[[cc]][[pp]][[kk]])),
                                               function(x) quantile(x, probs = seq(0, 1, length = 1000))))
        }
        waspRan20[[cc]][[pp]] <- matrix(NA, nrow = nrow(tmp[[1]]), ncol = ncol(vbRan[[cc]][[pp]]))
        for (dd in 1:ncol(waspRan20[[cc]][[pp]])) {
            waspRan20[[cc]][[pp]][ , dd] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , dd])))
        }
        tmp <- list()
        for (kk in 1:20) {
            tmp[[kk]] <- do.call(cbind, lapply(split(waspFixSamp[[cc]][[pp]][[kk]], col(waspFixSamp[[cc]][[pp]][[kk]])),
                                               function(x) quantile(x, probs = seq(0, 1, length = 1000))))
        }
        waspFix20[[cc]][[pp]] <- matrix(NA, nrow = nrow(tmp[[1]]), ncol = ncol(vbFix[[cc]][[pp]]))
        for (dd in 1:ncol(waspFix20[[cc]][[pp]])) {
            waspFix20[[cc]][[pp]][ , dd] <- rowMeans(do.call(cbind, lapply(tmp, function(x) x[ , dd])))
        }
    }
}

waspRan <- list()
waspFix <- list()
for (cc in 1:10) {
    waspRan[[cc]] <- list()
    waspFix[[cc]] <- list()
    for (pp in 1:2) {
        waspRan[[cc]][[pp]] <- list()
        waspFix[[cc]][[pp]] <- list()
        waspRan[[cc]][[pp]][[1]] <- waspRan10[[cc]][[pp]]
        waspRan[[cc]][[pp]][[2]] <- waspRan20[[cc]][[pp]]
        waspFix[[cc]][[pp]][[1]] <- waspFix10[[cc]][[pp]]
        waspFix[[cc]][[pp]][[2]] <- waspFix20[[cc]][[pp]]
    }
}

for (cc in 1:10) {
    waspTime[[cc]] <- list()
    for (pp in 1:2) {
        waspTime[[cc]][[pp]] <- list()
        waspTime[[cc]][[pp]][[1]] <- mean(unlist(waspTime10[[cc]][[pp]]))
        waspTime[[cc]][[pp]][[2]] <- mean(unlist(waspTime20[[cc]][[pp]]))
    }
}

for (cc in 1:10) {
    xingFix[[cc]] <- list()
    xingRan[[cc]] <- list()
    xingTime[[cc]] <- list()
    for (pp in 1:2) {
        xingFix[[cc]][[pp]] <- list()
        xingRan[[cc]][[pp]] <- list()
        xingTime[[cc]][[pp]] <- list()
        for (kk in 1:2) {
            if (kk == 1) {
                dat <- readRDS(paste0("xing/marg/xing_fix_ran_cv_", cc, "_p_", ndim[pp], "_k10.rds"))
            } else {
                dat <- readRDS(paste0("xing/marg/xing_fix_ran_cv_", cc, "_p_", ndim[pp], "_k20.rds"))
            }
            xingFix[[cc]][[pp]][[kk]] <- dat$fix
            xingRan[[cc]][[pp]][[kk]] <- dat$ran
            xingTime[[cc]][[pp]][[kk]] <- dat$time
        }
    }
}

for (cc in 1:10) {
    scotFix[[cc]] <- list()
    scotRan[[cc]] <- list()
    scotTime[[cc]] <- list()
    for (pp in 1:2) {
        scotFix[[cc]][[pp]] <- list()
        scotRan[[cc]][[pp]] <- list()
        scotTime[[cc]][[pp]] <- list()
        for (kk in 1:2) {
            if (kk == 1) {
                dat <- readRDS(paste0("cons/marg/cons_fix_ran_cv_", cc, "_p_", ndim[pp], "_k10.rds"))
            } else {
                dat <- readRDS(paste0("cons/marg/cons_fix_ran_cv_", cc, "_p_", ndim[pp], "_k20.rds"))
            }
            scotFix[[cc]][[pp]][[kk]] <- dat$fix
            scotRan[[cc]][[pp]][[kk]] <- dat$ran
            scotTime[[cc]][[pp]][[kk]] <- dat$time
        }
    }
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
    for (pp in 1:2) {
        resRan[["MCMC"]][[cc]][[pp]] <- list()
        resRan[["Xing"]][[cc]][[pp]] <- list()
        resRan[["Scot"]][[cc]][[pp]] <- list()
        resRan[["WASP"]][[cc]][[pp]] <- list()
        resRan[["VB"]][[cc]][[pp]] <- list()
        resRan[["bbvb"]][[cc]][[pp]] <- list()
        for (kk in 1:2) {
            resRan[["MCMC"]][[cc]][[pp]][[kk]] <- list()
            resRan[["Xing"]][[cc]][[pp]][[kk]] <- list()
            resRan[["Scot"]][[cc]][[pp]][[kk]] <- list()
            resRan[["WASP"]][[cc]][[pp]][[kk]] <- list()
            resRan[["VB"]][[cc]][[pp]][[kk]] <- list()
            resRan[["bbvb"]][[cc]][[pp]][[kk]] <- list()
            for (dd in 1:ncol(mcmcRan[[cc]][[pp]])) {
                rr <- range(c(mcmcRan[[cc]][[pp]][ , dd],
                              xingRan[[cc]][[pp]][[kk]][ , dd],
                              scotRan[[cc]][[pp]][[kk]][ , dd],
                              waspRan[[cc]][[pp]][[kk]][ , dd],
                              vbRan[[cc]][[pp]][ , dd],
                              bbvbRan[[cc]][[pp]][ , dd]
                              ))
                bw1 <- dpik(mcmcRan[[cc]][[pp]][ , dd], range.x = rr)
                bw2 <- dpik(xingRan[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw3 <- dpik(scotRan[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw4 <- dpik(waspRan[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw5 <- dpik(vbRan[[cc]][[pp]][ , dd], range.x = rr)
                bw6 <- dpik(bbvbRan[[cc]][[pp]][ , dd], range.x = rr)
                dens1 <- bkde(mcmcRan[[cc]][[pp]][ , dd], bandwidth = bw1, range.x = rr)
                dens2 <- bkde(xingRan[[cc]][[pp]][[kk]][ , dd], bandwidth = bw2, range.x = rr)
                dens3 <- bkde(scotRan[[cc]][[pp]][[kk]][ , dd], bandwidth = bw3, range.x = rr)
                dens4 <- bkde(waspRan[[cc]][[pp]][[kk]][ , dd], bandwidth = bw4, range.x = rr)
                dens5 <- bkde(vbRan[[cc]][[pp]][ , dd], bandwidth = bw5, range.x = rr)
                dens6 <- bkde(bbvbRan[[cc]][[pp]][ , dd], bandwidth = bw6, range.x = rr)
                resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]] <- dens1
                resRan[["Xing"]][[cc]][[pp]][[kk]][[dd]] <- dens2
                resRan[["Scot"]][[cc]][[pp]][[kk]][[dd]] <- dens3
                resRan[["WASP"]][[cc]][[pp]][[kk]][[dd]] <- dens4
                resRan[["VB"]][[cc]][[pp]][[kk]][[dd]] <- dens5
                resRan[["bbvb"]][[cc]][[pp]][[kk]][[dd]] <- dens6
            }
        }
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
    for (pp in 1:2) {
        resFix[["MCMC"]][[cc]][[pp]] <- list()
        resFix[["Xing"]][[cc]][[pp]] <- list()
        resFix[["Scot"]][[cc]][[pp]] <- list()
        resFix[["WASP"]][[cc]][[pp]] <- list()
        resFix[["VB"]][[cc]][[pp]] <- list()
        resFix[["bbvb"]][[cc]][[pp]] <- list()
        for (kk in 1:2) {
            resFix[["MCMC"]][[cc]][[pp]][[kk]] <- list()
            resFix[["Xing"]][[cc]][[pp]][[kk]] <- list()
            resFix[["Scot"]][[cc]][[pp]][[kk]] <- list()
            resFix[["WASP"]][[cc]][[pp]][[kk]] <- list()
            resFix[["VB"]][[cc]][[pp]][[kk]] <- list()
            resFix[["bbvb"]][[cc]][[pp]][[kk]] <- list()
            for (dd in 1:ncol(mcmcFix[[cc]][[pp]])) {
                rr <- range(c(mcmcFix[[cc]][[pp]][ , dd],
                              xingFix[[cc]][[pp]][[kk]][ , dd],
                              scotFix[[cc]][[pp]][[kk]][ , dd],
                              waspFix[[cc]][[pp]][[kk]][ , dd],
                              vbFix[[cc]][[pp]][ , dd],
                              bbvbFix[[cc]][[pp]][ , dd]
                              ))
                bw1 <- dpik(mcmcFix[[cc]][[pp]][ , dd], range.x = rr)
                bw2 <- dpik(xingFix[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw3 <- dpik(scotFix[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw4 <- dpik(waspFix[[cc]][[pp]][[kk]][ , dd], range.x = rr)
                bw5 <- dpik(vbFix[[cc]][[pp]][ , dd], range.x = rr)
                bw6 <- dpik(bbvbFix[[cc]][[pp]][ , dd], range.x = rr)
                dens1 <- bkde(mcmcFix[[cc]][[pp]][ , dd], bandwidth = bw1, range.x = rr)
                dens2 <- bkde(xingFix[[cc]][[pp]][[kk]][ , dd], bandwidth = bw2, range.x = rr)
                dens3 <- bkde(scotFix[[cc]][[pp]][[kk]][ , dd], bandwidth = bw3, range.x = rr)
                dens4 <- bkde(waspFix[[cc]][[pp]][[kk]][ , dd], bandwidth = bw4, range.x = rr)
                dens5 <- bkde(vbFix[[cc]][[pp]][ , dd], bandwidth = bw5, range.x = rr)
                dens6 <- bkde(bbvbFix[[cc]][[pp]][ , dd], bandwidth = bw6, range.x = rr)
                resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]] <- dens1
                resFix[["Xing"]][[cc]][[pp]][[kk]][[dd]] <- dens2
                resFix[["Scot"]][[cc]][[pp]][[kk]][[dd]] <- dens3
                resFix[["WASP"]][[cc]][[pp]][[kk]][[dd]] <- dens4
                resFix[["VB"]][[cc]][[pp]][[kk]][[dd]] <- dens5
                resFix[["bbvb"]][[cc]][[pp]][[kk]][[dd]] <- dens6
            }
        }
    }
}

accRan <- list()
for (pp in 1:2) {
    accRan[[pp]] <- list()
    for (kk in 1:2) {
        accRan[[pp]][[kk]] <- array(NA,
                                    dim = c(10, ncol = ncol(mcmcRan[[1]][[pp]]), 5),
                                    dimnames = list(paste0("cv", 1:10),
                                                    paste0("dim", 1:ncol(mcmcRan[[1]][[pp]])),
                                                    c("xing", "scot", "wasp", "vb", "bbvb")
                                                    )
                                    )
        for (cc in 1:10) {
            for (dd in 1:ncol(mcmcRan[[cc]][[pp]])) {
                accRan[[pp]][[kk]][cc, dd, 1] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resRan[["Xing"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accRan[[pp]][[kk]][cc, dd, 2] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resRan[["Scot"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accRan[[pp]][[kk]][cc, dd, 3] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resRan[["WASP"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accRan[[pp]][[kk]][cc, dd, 4] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resRan[["VB"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accRan[[pp]][[kk]][cc, dd, 5] <- (1 - sum(abs(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resRan[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resRan[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
            }
        }
    }
}

accVar <- array("NA", c(2, 2, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("xing", "scot", "wasp", "vb", "bbvb")))
accCov <- array("NA", c(2, 2, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("xing", "scot", "wasp", "vb", "bbvb")))
arrRan <- array("NA", c(2, 2, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("xing", "scot", "wasp", "vb", "bbvb")))

for (pp in 1:2) {
    for (kk in 1:2) {
        if (pp == 1) {
            vidx <- c(1, 4, 6)
            accVar[kk, pp, 1] <- paste(format(round(mean(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), paste0("(", format(round(sd(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), ")"))
            accVar[kk, pp, 2] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 2]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 2]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 3] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 3]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 3]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 4] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 4]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 4]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 5] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 5]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 5]), 2), ")")), nsmall = 2)
        } else {
            vidx <- c(1, 7, 12, 16, 19, 21)
            accVar[kk, pp, 1] <- paste(format(round(mean(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), paste0("(", format(round(sd(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), ")"))
            accVar[kk, pp, 2] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 2]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 2]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 3] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 3]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 3]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 4] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 4]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 4]), 2), ")")), nsmall = 2)
            accVar[kk, pp, 5] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 5]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 5]), 2), ")")), nsmall = 2)
        }

    }
}

for (pp in 1:2) {
    for (kk in 1:2) {
        if (pp == 1) {
            vidx <- -c(1, 4, 6)
            accCov[kk, pp, 1] <- paste(format(round(mean(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), paste0("(", format(round(sd(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), ")"))
            accCov[kk, pp, 2] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 2]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 2]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 3] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 3]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 3]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 4] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 4]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 4]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 5] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 5]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 5]), 2), ")")), nsmall = 2)
        } else {
            vidx <- -c(1, 7, 12, 16, 19, 21)
            accCov[kk, pp, 1] <- paste(format(round(mean(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), paste0("(", format(round(sd(accRan[[pp]][[kk]][ , vidx, 1]), 2), nsmall = 2), ")"))
            accCov[kk, pp, 2] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 2]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 2]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 3] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 3]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 3]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 4] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 4]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 4]), 2), ")")), nsmall = 2)
            accCov[kk, pp, 5] <- format(paste(round(mean(accRan[[pp]][[kk]][ , vidx, 5]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , vidx, 5]), 2), ")")), nsmall = 2)
        }

    }
}

arrRan <- array("NA", c(2, 2, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("xing", "scot", "wasp", "vb", "bbvb")))
for (pp in 1:2) {
    for (kk in 1:2) {
        arrRan[kk, pp, 1] <- format(paste(round(mean(accRan[[pp]][[kk]][ , , 1]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , , 1]), 2), ")")), nsmall = 2)
        arrRan[kk, pp, 2] <- format(paste(round(mean(accRan[[pp]][[kk]][ , , 2]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , , 2]), 2), ")")), nsmall = 2)
        arrRan[kk, pp, 3] <- format(paste(round(mean(accRan[[pp]][[kk]][ , , 3]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , , 3]), 2), ")")), nsmall = 2)
        arrRan[kk, pp, 4] <- format(paste(round(mean(accRan[[pp]][[kk]][ , , 4]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , , 4]), 2), ")")), nsmall = 2)
        arrRan[kk, pp, 5] <- format(paste(round(mean(accRan[[pp]][[kk]][ , , 5]), 2), paste0("(", round(sd(accRan[[pp]][[kk]][ , , 5]), 2), ")")), nsmall = 2)
    }
}

accFix <- list()
for (pp in 1:2) {
    accFix[[pp]] <- list()
    for (kk in 1:2) {
        accFix[[pp]][[kk]] <- array(NA,
                                    dim = c(10, ncol = ncol(mcmcFix[[1]][[pp]]), 5),
                                    dimnames = list(paste0("cv", 1:10),
                                                    paste0("dim", 1:ncol(mcmcFix[[1]][[pp]])),
                                                    c("xing", "scot", "wasp", "vb", "bbvb")
                                                    )
                                    )
        for (cc in 1:10) {
            for (dd in 1:ncol(mcmcFix[[cc]][[pp]])) {
                accFix[[pp]][[kk]][cc, dd, 1] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resFix[["Xing"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accFix[[pp]][[kk]][cc, dd, 2] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resFix[["Scot"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accFix[[pp]][[kk]][cc, dd, 3] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resFix[["WASP"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accFix[[pp]][[kk]][cc, dd, 4] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resFix[["VB"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
                accFix[[pp]][[kk]][cc, dd, 5] <- (1 - sum(abs(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$y  - resFix[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$y) * diff(resFix[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x)[1]) / 2)
            }
        }
    }
}

arrFix <- array("NA", c(2, 2, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("xing", "scot", "wasp", "vb", "bbvb")))
for (pp in 1:2) {
    for (kk in 1:2) {
        arrFix[kk, pp, 1] <- format(paste(round(mean(accFix[[pp]][[kk]][ , , 1]), 2), paste0("(", round(sd(accFix[[pp]][[kk]][ , , 1]), 2), ")")), nsmall = 2)
        arrFix[kk, pp, 2] <- format(paste(round(mean(accFix[[pp]][[kk]][ , , 2]), 2), paste0("(", round(sd(accFix[[pp]][[kk]][ , , 2]), 2), ")")), nsmall = 2)
        arrFix[kk, pp, 3] <- format(paste(round(mean(accFix[[pp]][[kk]][ , , 3]), 2), paste0("(", round(sd(accFix[[pp]][[kk]][ , , 3]), 2), ")")), nsmall = 2)
        arrFix[kk, pp, 4] <- format(paste(round(mean(accFix[[pp]][[kk]][ , , 4]), 2), paste0("(", round(sd(accFix[[pp]][[kk]][ , , 4]), 2), ")")), nsmall = 2)
        arrFix[kk, pp, 5] <- format(paste(round(mean(accFix[[pp]][[kk]][ , , 5]), 2), paste0("(", round(sd(accFix[[pp]][[kk]][ , , 5]), 2), ")")), nsmall = 2)
    }
}

fixTbl <- rbind(as.vector(arrFix[ , , "scot"]), as.vector(arrFix[ , , "xing"]), as.vector(arrFix[ , , "vb"]),
                as.vector(arrFix[ , , "bbvb"]), as.vector(arrFix[ , , "wasp"]))
ranTbl <- rbind(as.vector(arrRan[ , , "scot"]), as.vector(arrRan[ , , "xing"]), as.vector(arrRan[ , , "vb"]),
                as.vector(arrRan[ , , "bbvb"]), as.vector(arrRan[ , , "wasp"]))
rownames(fixTbl) <- rownames(ranTbl) <- c("CMC", "SDP", "VB", "BBVB", "WASP")
colnames(fixTbl) <- colnames(ranTbl) <- c("p10, k10", "p10, k20", "p100, k10", "p100, k20")

varTbl <- rbind(as.vector(accVar[ , , "scot"]), as.vector(accVar[ , , "xing"]),
                as.vector(accVar[ , , "vb"]), as.vector(accVar[ , , "bbvb"]),
                as.vector(accVar[ , , "wasp"]))
covTbl <- rbind(as.vector(accCov[ , , "scot"]), as.vector(accCov[ , , "xing"]), as.vector(accCov[ , , "vb"]),
                as.vector(accCov[ , , "bbvb"]), as.vector(accCov[ , , "wasp"]))
rownames(varTbl) <- rownames(covTbl) <- c("CMC", "SDP", "VB", "BBVB", "WASP")
colnames(varTbl) <- colnames(covTbl) <- c("p10, k10", "p10, k20", "p100, k10", "p100, k20")

xtable::xtable(fixTbl)
xtable::xtable(varTbl)
xtable::xtable(covTbl)

### 2d plots
mcmcCovMat <- list()
vbCovMat <- list()
bbvbCovMat <- list()
scotCovMat <- list()
xingCovMat <- list()
waspCovMat <- list()

for (cc in 1:10) {
    mcmcCovMat[[cc]] <- list()
    vbCovMat[[cc]] <- list()
    bbvbCovMat[[cc]] <- list()
    for (pp in 1:2) {
        mcmcCovMat[[cc]][[pp]] <- list()
        vbCovMat[[cc]][[pp]] <- list()
        bbvbCovMat[[cc]][[pp]] <- list()
        for (dd in 1:3) {
            if (pp == 1) {
                cov2d <- cbind(c(2, 2, 3), c(3, 5, 5))
                mcmcCovMat[[cc]][[pp]][[dd]] <- mcmcRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                vbCovMat[[cc]][[pp]][[dd]] <- vbRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                bbvbCovMat[[cc]][[pp]][[dd]] <- bbvbRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
            } else {
                cov2d <- cbind(c(2, 2, 3), c(3, 8, 8))
                mcmcCovMat[[cc]][[pp]][[dd]] <- mcmcRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                vbCovMat[[cc]][[pp]][[dd]] <- vbRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                bbvbCovMat[[cc]][[pp]][[dd]] <- bbvbRan[[cc]][[pp]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
            }
        }
    }
}

ndim <- 1:2
for (cc in 1:10) {
    xingCovMat[[cc]] <- list()
    for (pp in 1:2) {
        xingCovMat[[cc]][[pp]] <- list()
        for (kk in 1:2) {
            if (kk == 1) {
                cdat <- readRDS(paste0("xing/joint/xing_cov_cv_", cc, "_p_", ndim[pp], "_k10.rds"))
            } else {
                cdat <- readRDS(paste0("xing/joint/xing_cov_cv_", cc, "_p_", ndim[pp], "_k20.rds"))
            }
            xingCovMat[[cc]][[pp]][[kk]] <- list()
            for (dd in 1:3) {
                xingCovMat[[cc]][[pp]][[kk]][[dd]] <- cdat$cov[[dd]]
            }
        }
    }
}

ndim <- 1:2
for (cc in 1:10) {
    scotCovMat[[cc]] <- list()
    for (pp in 1:2) {
        scotCovMat[[cc]][[pp]] <- list()
        for (kk in 1:2) {
            if (kk == 1) {
                cdat <- readRDS(paste0("cons/joint/cons_cov_cv_", cc, "_p_", ndim[pp], "_k10.rds"))
            } else {
                cdat <- readRDS(paste0("cons/joint/cons_cov_cv_", cc, "_p_", ndim[pp], "_k20.rds"))
            }
            scotCovMat[[cc]][[pp]][[kk]] <- list()
            for (dd in 1:3) {
                scotCovMat[[cc]][[pp]][[kk]][[dd]] <- cdat$cov[[dd]]
            }
        }
    }
}

for (cc in 1:10) {
    waspCovMat[[cc]] <- list()
    for (pp in 1:2) {
        waspCovMat[[cc]][[pp]] <- list()
        for (kk in 1:2) {
            waspCovMat[[cc]][[pp]][[kk]] <- list()
            for (dd in 1:3) {
                if (kk == 1) {
                    dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_p_", pp, "_d_", dd, "_k10.csv"), header = FALSE)
                    pr <- as.numeric(dat[ , 3])
                    pr[pr < 1e-10] <- 0.0
                    waspCovMat[[cc]][[pp]][[kk]][[dd]] <- dat[sample(1:nrow(dat), 10000, replace = TRUE, prob = pr), 1:2]
                } else {
                    dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_p_", pp, "_d_", dd, "_k20.csv"), header = FALSE)
                    pr <- as.numeric(dat[ , 3])
                    pr[pr < 1e-10] <- 0.0
                    waspCovMat[[cc]][[pp]][[kk]][[dd]] <- dat[sample(1:nrow(dat), 10000, replace = TRUE, prob = pr), 1:2]
                }

            }
        }
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
    for (pp in 1:2) {
        resCovMat[["MCMC"]][[cc]][[pp]] <- list()
        resCovMat[["Xing"]][[cc]][[pp]] <- list()
        resCovMat[["Scot"]][[cc]][[pp]] <- list()
        resCovMat[["WASP"]][[cc]][[pp]] <- list()
        resCovMat[["VB"]][[cc]][[pp]] <- list()
        resCovMat[["bbvb"]][[cc]][[pp]] <- list()
        for (kk in 1:2) {
            resCovMat[["MCMC"]][[cc]][[pp]][[kk]] <- list()
            resCovMat[["Xing"]][[cc]][[pp]][[kk]] <- list()
            resCovMat[["Scot"]][[cc]][[pp]][[kk]] <- list()
            resCovMat[["WASP"]][[cc]][[pp]][[kk]] <- list()
            resCovMat[["VB"]][[cc]][[pp]][[kk]] <- list()
            resCovMat[["bbvb"]][[cc]][[pp]][[kk]] <- list()
            for (dd in 1:3) {
                if (kk == 1) {
                    dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_p_", pp, "_d_", dd, "_k20.csv"), header = FALSE)
                    bw1 <- max(diff(as.numeric(dat[ , 1])))
                    bw2 <- max(diff(as.numeric(dat[ , 2])))
                } else {
                    dat <- read.csv(paste0("wasp/joint/wasp_cov_cv_", cc, "_p_", pp, "_d_", dd, "_k20.csv"), header = FALSE)
                    bw1 <- max(diff(as.numeric(dat[ , 1])))
                    bw2 <- max(diff(as.numeric(dat[ , 2])))
                }
                rr1 <- range(c(mcmcCovMat[[cc]][[pp]][[dd]][ , 1],
                              xingCovMat[[cc]][[pp]][[kk]][[dd]][ , 1],
                              scotCovMat[[cc]][[pp]][[kk]][[dd]][ , 1],
                              waspCovMat[[cc]][[pp]][[kk]][[dd]][ , 1],
                              vbCovMat[[cc]][[pp]][[dd]][ , 1],
                              bbvbCovMat[[cc]][[pp]][[dd]][ , 1]
                              ))
                rr2 <- range(c(mcmcCovMat[[cc]][[pp]][[dd]][ , 2],
                              xingCovMat[[cc]][[pp]][[kk]][[dd]][ , 2],
                              scotCovMat[[cc]][[pp]][[kk]][[dd]][ , 2],
                              waspCovMat[[cc]][[pp]][[kk]][[dd]][ , 2],
                              vbCovMat[[cc]][[pp]][[dd]][ , 2],
                              bbvbCovMat[[cc]][[pp]][[dd]][ , 2]
                              ))
                bw12 <- bw22 <- bw32 <- bw42 <- bw52 <- bw62 <- bw2
                bw11 <- bw21 <- bw31 <- bw41 <- bw51 <- bw61 <- bw1
                dens1 <- bkde2D(mcmcCovMat[[cc]][[pp]][[dd]], bandwidth = c(bw11, bw12), range.x = list(rr1, rr2))
                dens2 <- bkde2D(xingCovMat[[cc]][[pp]][[kk]][[dd]], bandwidth = c(bw21, bw22), range.x = list(rr1, rr2))
                dens3 <- bkde2D(scotCovMat[[cc]][[pp]][[kk]][[dd]], bandwidth = c(bw31, bw32), range.x = list(rr1, rr2))
                dens4 <- bkde2D(waspCovMat[[cc]][[pp]][[kk]][[dd]], bandwidth = c(bw41, bw42), range.x = list(rr1, rr2))
                dens5 <- bkde2D(vbCovMat[[cc]][[pp]][[dd]], bandwidth = c(bw51, bw52), range.x = list(rr1, rr2))
                dens6 <- bkde2D(bbvbCovMat[[cc]][[pp]][[dd]], bandwidth = c(bw61, bw62), range.x = list(rr1, rr2))
                resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]] <- dens1
                resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]] <- dens2
                resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]] <- dens3
                resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]] <- dens4
                resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]] <- dens5
                resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]] <- dens6
            }
        }
    }
}

accCovMat <- list()
for (pp in 1:2) {
    accCovMat[[pp]] <- list()
    for (kk in 1:2) {
        accCovMat[[pp]][[kk]] <- array(NA,
                                    dim = c(10, 3, 5),
                                    dimnames = list(paste0("cv", 1:10),
                                                    paste0("dim", 1:3),
                                                    c("xing", "scot", "wasp", "vb", "bbvb")
                                                    )
                                    )
        for (cc in 1:10) {
            for (dd in 1:3) {
                accCovMat[[pp]][[kk]][cc, dd, 1] <- 1 - sum(abs(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat - resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2)[1]) / 2
                accCovMat[[pp]][[kk]][cc, dd, 2] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat - resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2)[1]) / 2
                accCovMat[[pp]][[kk]][cc, dd, 3] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat - resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2)[1]) / 2
                accCovMat[[pp]][[kk]][cc, dd, 4] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat - resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2)[1]) / 2
                accCovMat[[pp]][[kk]][cc, dd, 5] <-  1 - sum(abs(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat - resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$fhat) * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1)[1] * diff(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2)[1]) / 2
            }
        }
    }
}

arrCovMat <- array("NA", c(2, 2, 3, 5), dimnames = list(c("k10", "k20"), c("p10", "p100"), c("d1", "d2", "d3"), c("xing", "scot", "wasp", "vb", "bbvb")))
for (pp in 1:2) {
    for (dd in 1:3) {
        for (kk in 1:2) {
            arrCovMat[kk, pp, dd, 1] <- paste(format(round(mean(accCovMat[[pp]][[kk]][ , dd, 1]), 2), nsmall = 2),
                                               paste0("(", format(round(sd(accCovMat[[pp]][[kk]][ , dd, 1]), 2), nsmall = 2), ")"))
            arrCovMat[kk, pp, dd, 2] <- paste(format(round(mean(accCovMat[[pp]][[kk]][ , dd, 2]), 2), nsmall = 2),
                                               paste0("(", format(round(sd(accCovMat[[pp]][[kk]][ , dd, 2]), 2), nsmall = 2), ")"))
            arrCovMat[kk, pp, dd, 3] <- paste(format(round(mean(accCovMat[[pp]][[kk]][ , dd, 3]), 2), nsmall = 2),
                                               paste0("(", format(round(sd(accCovMat[[pp]][[kk]][ , dd, 3]), 2), nsmall = 2), ")"))
            arrCovMat[kk, pp, dd, 4] <- paste(format(round(mean(accCovMat[[pp]][[kk]][ , dd, 4]), 2), nsmall = 2),
                                              paste0("(", format(round(sd(accCovMat[[pp]][[kk]][ , dd, 4]), 2), nsmall = 2), ")"))
            arrCovMat[kk, pp, dd, 5] <- paste(format(round(mean(accCovMat[[pp]][[kk]][ , dd, 5]), 2), nsmall = 2),
                                               paste0("(", format(round(sd(accCovMat[[pp]][[kk]][ , dd, 5]), 2), nsmall = 2), ")"))
        }
    }
}

covMatTbl1 <- rbind(as.vector(arrCovMat[ , 1, , "scot"]), as.vector(arrCovMat[ , 1, , "xing"]),
                    as.vector(arrCovMat[ , 1, , "vb"]), as.vector(arrCovMat[ , 1, , "bbvb"]),
                    as.vector(arrCovMat[ , 1, , "wasp"]))
rownames(covMatTbl1) <- c("CMC", "SDP", "VB", "BBVB", "WASP")
colnames(covMatTbl1) <-  paste(rep(paste0("d", 1:3), each = 2), rep(c("k10", "k20"), times = 3), sep = ", ")

covMatTbl2 <- rbind(as.vector(arrCovMat[ , 2, , "scot"]), as.vector(arrCovMat[ , 2, , "xing"]),
                    as.vector(arrCovMat[ , 2, , "vb"]), as.vector(arrCovMat[ , 2, , "bbvb"]),
                    as.vector(arrCovMat[ , 2, , "wasp"]))
rownames(covMatTbl2) <- c("CMC", "SDP", "VB", "BBVB", "WASP")
colnames(covMatTbl2) <-  paste(rep(paste0("d", 1:3), each = 2), rep(c("k10", "k20"), times = 3), sep = ", ")

xtable::xtable(covMatTbl1)
xtable::xtable(covMatTbl2)

rnames <- list(c("12", "13"), c("12", "23"), c("13", "23"))

cc <- 1
pp <- 1

nparts <- c(10, 20)
pdf("~/wasp/lme/result/img/cov_2d_r3.pdf", width = 35, height = 15)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(2, 4, 0, 2), oma = c(1, 1, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (kk in 1:2) {
    for (dd in 1:3) {
        contour(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat,
                xlim = range(c(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x1,
                             resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x1)) + c(-0.01, 0.01),
                ylim = range(c(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x2,
                             resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x2)) + c(-0.01, 0.01),
                lwd = 6, axes = FALSE, col = colors[1]
                )
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        contour(resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[2])
        contour(resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[3])
        contour(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[4])
        contour(resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[5])
        contour(resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[6])
        axis(side = 2, tck = -0.01, lwd = 4, cex = 2, labels = NA)
        mtext(format(round(yylab, 3)), at = yylab, side = 2, line = 0.5, cex = 2.5, las = 2)
        axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
        grid(lwd = 2)
        box(col = "grey40", lwd = 4)
        mtext(paste0("k = ", bquote(.(nparts[kk])), ", r = 3"), side = 3, line = -7, adj = 0.1, cex = 3.5)
        mtext(bquote(sigma[.(rnames[[dd]][1])]), side = 3, line = -4, adj = 0.1, cex = 4)
        mtext(", ", side = 3, line = -3.5, adj = 0.22, cex = 4)
        mtext(bquote(sigma[.(rnames[[dd]][2])]), side = 3, line = -4, adj = 0.28, cex = 4)
        if (kk == 2 & dd == 2) {
            legend("bottom", c("ADVI", "CMC", "MCMC", "SA", "SDP", "WASP"), lty = c("solid"),
                   col = colors[c(6, 2, 1, 4, 3, 5)], lwd = 8, bty = "n", cex = 2.5, ncol = 3)
        }
    }
}
dev.off()

cc <- 1
pp <- 2
nparts <- c(10, 20)

pdf("~/wasp/lme/result/img/cov_2d_r6.pdf", width = 35, height = 15)
par(mfrow = c(2, 3))
par(cex = 1)
par(mar = c(2, 4, 0, 2), oma = c(1, 1, 0.6, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
for (kk in 1:2) {
    for (dd in 1:3) {
        contour(resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$fhat,
                xlim = range(c(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x1,
                               resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x1)),
                ylim = range(c(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["MCMC"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x2,
                               resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x2)),
                lwd = 6, axes = FALSE, col = colors[1]
                )
        xxlab <- axTicks(1)
        yylab <- axTicks(2)
        contour(resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Scot"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[2])
        contour(resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["Xing"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[3])
        contour(resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["VB"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[4])
        contour(resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["WASP"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[5])
        contour(resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x1, resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$x2, resCovMat[["bbvb"]][[cc]][[pp]][[kk]][[dd]]$fhat, add = TRUE, lwd = 6, col = colors[6])
        axis(side = 2, tck = -0.01, lwd = 4, cex = 2, labels = NA)
        mtext(format(round(yylab, 3)), at = yylab, side = 2, line = 0.5, cex = 2.5, las = 2)
        axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
        mtext(format(round(xxlab, 2)), at = round(xxlab, 2), side = 1, line = 1, cex = 2.5)
        grid(lwd = 2)
        box(col = "grey40", lwd = 4)
        mtext(paste0("k = ", bquote(.(nparts[kk])), ", r = 6"), side = 3, line = -7, adj = 0.1, cex = 3.5)
        mtext(bquote(sigma[.(rnames[[dd]][1])]), side = 3, line = -4, adj = 0.1, cex = 4)
        mtext(", ", side = 3, line = -3.5, adj = 0.22, cex = 4)
        mtext(bquote(sigma[.(rnames[[dd]][2])]), side = 3, line = -4, adj = 0.28, cex = 4)

        if (kk == 2 & dd == 2) {
            legend("bottom", c("ADVI", "CMC", "MCMC", "SA", "SDP", "WASP"), lty = c("solid"),
                   col = colors[c(6, 2, 1, 4, 3, 5)], lwd = 8, bty = "n", cex = 2.5, ncol = 3)
        }
    }
}
dev.off()

nparts <- c(10, 20)
for (cc in 1:10) {
    for (pp in 1:2) {
        for (kk in 1:2) {
            dat <- read.csv(paste0("wasp/joint/cov_2d_times_cv_", cc, "_k", nparts[kk], ".csv"), header = FALSE)
            waspTime[[cc]][[pp]][[kk]] <- waspTime[[cc]][[pp]][[kk]] + mean(unlist(dat[ , pp]))
        }
    }
}

rtime100 <- list(log10(unlist(lapply(mcmcTime, function(x) x[[1]]))),
                 log10(unlist(lapply(bbvbTime, function(x) x[[1]]))),
                 log10(unlist(lapply(vbTime, function(x) x[[1]]))),
                 log10(unlist(lapply(scotTime, function(x) x[[1]][[1]]))),
                 log10(unlist(lapply(scotTime, function(x) x[[1]][[2]]))),
                 log10(unlist(lapply(xingTime, function(x) x[[1]][[1]]))),
                 log10(unlist(lapply(xingTime, function(x) x[[1]][[2]]))),
                 log10(unlist(lapply(waspTime, function(x) x[[1]][[1]]))),
                 log10(unlist(lapply(waspTime, function(x) x[[1]][[2]])))
                 )

rtime200 <- list(log10(unlist(lapply(mcmcTime, function(x) x[[2]]))),
                 log10(unlist(lapply(bbvbTime, function(x) x[[2]]))),
                 log10(unlist(lapply(vbTime, function(x) x[[2]]))),
                 log10(unlist(lapply(scotTime, function(x) x[[2]][[1]]))),
                 log10(unlist(lapply(scotTime, function(x) x[[2]][[2]]))),
                 log10(unlist(lapply(xingTime, function(x) x[[2]][[1]]))),
                 log10(unlist(lapply(xingTime, function(x) x[[2]][[2]]))),
                 log10(unlist(lapply(waspTime, function(x) x[[2]][[1]]))),
                 log10(unlist(lapply(waspTime, function(x) x[[2]][[2]])))
                 )


pdf("~/wasp/lme/result/img/lme_time.pdf", 10, 8)
par(mfrow = c(1, 2))
par(cex = 1)
par(mar = c(6.7, 0, 0, 0), oma = c(3, 5, 0.2, 0.2))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))
for (pp in 1:2) {
    if (pp == 1) {
        boxplot(rtime100, ylab = NA, axes = FALSE, lwd = 2, ylim = c(1, 6),
                boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
                staplelwd = 3, medlwd = 4, outcex = 1.5)
        text(labels = c("p = 4, q = 3"), x = 5, y = 5.75, cex = 3)
    } else {
        boxplot(rtime200, ylab = NA, axes = FALSE, lwd = 2, ylim = c(1, 6),
                boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
                staplelwd = 3, medlwd = 4, outcex = 1.5)
        text(labels = c("p = 80, q = 6"), x = 5, y = 5.75, cex = 3)
    }
    grid(lwd=3)
    box(col = "grey40", lwd = 4)
    axis(side = 1, tck = -.01, labels = NA)
    mtext(c("MCMC", "ADVI", "SA", "CMC (k=10)", "CMC (k=20)", "SDP (k=10)",
            "SDP (k=20)", "WASP (k=10)", "WASP (k=20)"),
          at = c(1:9), side = 1, line = 0, cex = 1.9, las = 2)
    if (pp == 1) {
        axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
        mtext(format(seq(1, 6, by = 1)), at = seq(1, 6, by = 1), side = 2, line = 0.3, las = 1, cex = 2)
        mtext(expression(log[10] * " Seconds"), side = 2, outer = TRUE, cex = 2.5, line = 2.5)
    }
}
dev.off()
