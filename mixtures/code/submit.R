cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("full_sampler.R")
    cvtrain <- readRDS("../data/mix_100k.rds")
    train <- cvtrain[[id]]
    dataMat <- train$data
    res <- mvnMix(dataMat, ncomp = 2, niter = 10000, nburn = 5000, nthin = 5)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/full/res_", id, "_100k.rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("wasp_sampler.R")
    cvtrain <- readRDS("../data/wasp_mix_100k_k10.rds")
    reps <- rep(1:10, each = 10)
    subs <- rep(1:10, times = 10)

    wids <- cbind(reps, subs)
    cid <- wids[id, 1]
    sid <- wids[id, 2]

    train <- cvtrain[[cid]][[sid]]
    res <- mvnWaspMix(train, ncomp = 2, nrep = 10, niter = 10000, nburn = 5000, nthin = 5)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/sub10/samp/res_cv_", cid, "_nsub_", sid, "_k10.rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("wasp_sampler.R")
    cvtrain <- readRDS("../data/wasp_mix_100k_k5.rds")

    reps <- rep(1:10, each = 5)
    subs <- rep(1:5, times = 10)

    wids <- cbind(reps, subs)
    cid <- wids[id, 1]
    sid <- wids[id, 2]

    train <- cvtrain[[cid]][[sid]]

    res <- mvnWaspMix(train, ncomp = 2, nrep = 5, niter = 10000, nburn = 5000, nthin = 5)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/sub5/samp/res_cv_", cid, "_nsub_", sid, "_k5.rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("comp_sampler.R")
    cvtrain <- readRDS("../data/wasp_mix_100k_k10.rds")
    reps <- rep(1:10, each = 10)
    subs <- rep(1:10, times = 10)

    wids <- cbind(reps, subs)
    cid <- wids[id, 1]
    sid <- wids[id, 2]

    train <- cvtrain[[cid]][[sid]]
    res <- mvnCompMix(train, ncomp = 2, nrep = 10, niter = 10000, nburn = 5000, nthin = 5)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub10/comp_cv_", cid, "_nsub_", sid, "_k10.rds")
    saveRDS(res, fname)
} else if (mtd == 5) {
    source("comp_sampler.R")
    cvtrain <- readRDS("../data/wasp_mix_100k_k5.rds")

    reps <- rep(1:10, each = 5)
    subs <- rep(1:5, times = 10)

    wids <- cbind(reps, subs)
    cid <- wids[id, 1]
    sid <- wids[id, 2]

    train <- cvtrain[[cid]][[sid]]

    res <- mvnCompMix(train, ncomp = 2, nrep = 5, niter = 10000, nburn = 5000, nthin = 5)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub5/comp_cv_", cid, "_nsub_", sid, "_k5.rds")
    saveRDS(res, fname)
} else if (mtd == 6) {
    source("vb_sampler.R")

    cvtrain <- readRDS("../data/mix_100k.rds")
    train <- cvtrain[[id]]
    dataMat <- train$data

    res <- mvnVbMix(dataMat, ncomp = 2, niter = 1000)
    fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/vb/vb_cv_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 7) {
    library(parallelMCMCcombine)

    cvs <- rep(1:10, each = 2)
    subs <- rep(1:2, times = 10)
    cid <- cvs[id]
    sid <- subs[id]

    xx <- seq(0, 10, length = 500)
    if (sid == 1) {
        subdens <- array(0.0, dim = c(500, 1000, 5))
        tmp <- numeric(5)
        for (kk in 1:5) {
            fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub5/comp_cv_", cid, "_nsub_", kk, "_k5.rds")
            samp <- readRDS(fname)
            cat("kk: ", kk, "\n")
            for (ii in 1:1000) {
                for (gg in seq_along(xx)) {
                    yy <- c(xx[gg], xx[gg])
                    dens1 <- dmvnorm(yy, mean = samp$mu[1, , ii], sigma = samp$cov[1, , , ii])
                    dens2 <- dmvnorm(yy, mean = samp$mu[2, , ii], sigma = samp$cov[2, , , ii])
                    subdens[gg, ii, kk] <- samp$prob[ii, 1] * dens1 + samp$prob[ii, 2] * dens2
                }
            }
            tmp[kk] <- samp$time
        }
        fname1 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/cons/cons_dens_cv_", cid, "_k5.rds")
        fname2 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/xing/xing_dens_cv_", cid, "_k5.rds")
    } else {
        subdens <- array(0.0, dim = c(500, 1000, 10))
        tmp <- numeric(10)
        for (kk in 1:10) {
            fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub10/comp_cv_", cid, "_nsub_", kk, "_k10.rds")
            samp <- readRDS(fname)
            for (ii in 1:1000) {
                for (gg in seq_along(xx)) {
                    yy <- c(xx[gg], xx[gg])
                    dens1 <- dmvnorm(yy, mean = samp$mu[1, , ii], sigma = samp$cov[1, , , ii])
                    dens2 <- dmvnorm(yy, mean = samp$mu[2, , ii], sigma = samp$cov[2, , , ii])
                    subdens[gg, ii, kk] <- samp$prob[ii, 1] * dens1 + samp$prob[ii, 2] * dens2
                }
            }
            tmp[kk] <- samp$time
        }
        fname1 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/cons/cons_dens_cv_", cid, "_k10.rds")
        fname2 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/xing/xing_dens_cv_", cid, "_k10.rds")
    }

    strt1 <- proc.time()
    scottDens <- consensusMCindep(subchain = subdens)
    end1 <- proc.time()
    strt2 <- proc.time()
    try1 <- tryCatch(xingDens <- semiparamDPE(subchain = subdens),
                     error = function(e) e)
    if (any(class(try1) == "simpleError")) {
        xingDens <- try1
    }
    end2 <- proc.time()

    saveRDS(list("dens" = t(scottDens), time = mean(tmp) + end1[3] - strt1[3]), fname1)
    saveRDS(list("dens" = t(xingDens), time = mean(tmp) + end2[3] - strt2[3]), fname2)
} else (mtd == 8) {
    library(parallelMCMCcombine)

    cvs <- rep(1:10, each = 2)
    subs <- rep(1:2, times = 10)
    cid <- cvs[id]
    sid <- subs[id]

    if (sid == 1) {
        subdensRho <- array(0.0, dim = c(2, 1000, 5))
        tmp <- numeric(5)
        for (kk in 1:5) {
            fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub5/comp_cv_", cid, "_nsub_", kk, "_k5.rds")
            samp <- readRDS(fname)
            ppp <- colMeans(samp$prob)
            if (ppp[1] < ppp[2]) {
                for (ii in 1:1000) {
                    subdensRho[1, ii, kk] <- cov2cor(samp$cov[1, , , ii])[1, 2]
                    subdensRho[2, ii, kk] <- cov2cor(samp$cov[2, , , ii])[1, 2]
                }
            } else {
                for (ii in 1:1000) {
                    subdensRho[2, ii, kk] <- cov2cor(samp$cov[1, , , ii])[1, 2]
                    subdensRho[1, ii, kk] <- cov2cor(samp$cov[2, , , ii])[1, 2]
                }
            }
            tmp[kk] <- samp$time
        }
        fname1 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/cons/cons_rho_cv_", cid, "_k5.rds")
        fname2 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/xing/xing_rho_cv_", cid, "_k5.rds")
    } else {
        subdensRho <- array(0.0, dim = c(2, 1000, 10))
        tmp <- numeric(10)
        for (kk in 1:10) {
            fname <- paste0("/Shared/ssrivastva/wasp/mixtures/result/comp/sub10/comp_cv_", cid, "_nsub_", kk, "_k10.rds")
            samp <- readRDS(fname)
            ppp <- colMeans(samp$prob)
            if (ppp[1] < ppp[2]) {
                for (ii in 1:1000) {
                    subdensRho[1, ii, kk] <- cov2cor(samp$cov[1, , , ii])[1, 2]
                    subdensRho[2, ii, kk] <- cov2cor(samp$cov[2, , , ii])[1, 2]
                }
            } else {
                for (ii in 1:1000) {
                    subdensRho[2, ii, kk] <- cov2cor(samp$cov[1, , , ii])[1, 2]
                    subdensRho[1, ii, kk] <- cov2cor(samp$cov[2, , , ii])[1, 2]
                }
            }
            tmp[kk] <- samp$time
        }
        fname1 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/cons/cons_rho_cv_", cid, "_k10.rds")
        fname2 <- paste0("/Shared/ssrivastva/wasp/mixtures/result/xing/xing_rho_cv_", cid, "_k10.rds")
    }

    strt1 <- proc.time()
    scottDens <- consensusMCindep(subchain = subdensRho)
    end1 <- proc.time()
    strt2 <- proc.time()
    try1 <- tryCatch(xingDens <- semiparamDPE(subchain = subdensRho),
                     error = function(e) e)
    if (any(class(try1) == "simpleError")) {
        xingDens <- try1
    }
    end2 <- proc.time()

    saveRDS(list("dens" = t(scottDens), time = mean(tmp) + end1[3] - strt1[3]), fname1)
    saveRDS(list("dens" = t(xingDens), time = mean(tmp) + end2[3] - strt2[3]), fname2)
}
