cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("mcmc_sampler.R")
    cvtrain <- readRDS("../data/ml_train.rds")
    train <- cvtrain[[id]]
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)
    res <- sampleFromMixMdl(yvec, xmat, zmat, group, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/ml/result/full/full_res_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("variational_bayes.R")
    cvtrain <- readRDS("../data/ml_train.rds")
    train <- cvtrain[[id]]
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)
    res <- fitLinearMixefEffectsVB(yvec, xmat, zmat, group, 1000)
    fname <- paste0("/Shared/ssrivastva/wasp/ml/result/vb/vb_res_", id, ".rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("wasp_sampler.R")
    cvs <- rep(1:10, each = 10)
    wids <- cbind(cvs, rep(1:10, times = 10))

    cid <- wids[id, 1]
    sid <- wids[id, 2]

    cvtrain <- readRDS("../data/wasp_ml_train.rds")
    train <- cvtrain[[cid]][[sid]]
    rm(cvtrain)
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)
    res <- sampleFromWaspMixMdl(yvec, xmat, zmat, group, as.numeric(train$nrep), 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/ml/result/wasp/samp/wasp_cv_", cid, "_sub_", sid, "_k10.rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("comp_sampler.R")
    cvs <- rep(1:10, each = 10)
    wids <- cbind(cvs, rep(1:10, times = 10))

    cid <- wids[id, 1]
    sid <- wids[id, 2]

    cvtrain <- readRDS("../data/wasp_ml_train.rds")
    train <- cvtrain[[cid]][[sid]]
    rm(cvtrain)
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)
    res <- sampleFromCompMixMdl (yvec, xmat, zmat, group, as.numeric(train$nrep), 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/ml/result/comp/samp/comp_cv_", cid, "_sub_", sid, "_k10.rds")
    saveRDS(res, fname)
} else if (mtd == 5) {
    library(parallelMCMCcombine)

    cid <- id

    subfix <- array(0.0, dim = c(6, 1000, 10))
    subran <- array(0.0, dim = c(21, 1000, 10))
    tmp <- numeric(10)

    meanMap <- 1:6
    covMap <- c(7:12, 14:18, 21:24, 28:30, 35:36, 42)
    for (kk in 1:10) {
        fname <- paste0("/Shared/ssrivastva/wasp/ml/result/comp/samp/comp_cv_", cid, "_sub_", kk, "_k10.rds")
        samp <- readRDS(fname)
        subfix[ , , kk] <- t(samp$samples[ , meanMap])
        subran[ , , kk] <- t(samp$samples[ , covMap])
        tmp[kk] <- samp$time[3]
    }

    stime <- rep(0, 2)
    strt1 <- proc.time()
    scottFix <- consensusMCindep(subchain = subfix)
    end1 <- proc.time()
    stime[1] <- mean(tmp) + end1[3] - strt1[3]
    strt1 <- proc.time()
    scottRan <- consensusMCindep(subchain = subran)
    end1 <- proc.time()
    stime[2] <- mean(tmp) + end1[3] - strt1[3]

    xtime <- rep(0, 2)
    strt2 <- proc.time()
    xingFix <- semiparamDPE(subchain = subfix)
    end2 <- proc.time()
    xtime[1] <- mean(tmp) + end2[3] - strt2[3]
    strt2 <- proc.time()
    xingRan <- semiparamDPE(subchain = subran)
    end2 <- proc.time()
    xtime[2] <- mean(tmp) + end2[3] - strt2[3]

    fname1 <- paste0("/Shared/ssrivastva/wasp/ml/result/cons/marg/cons_fix_ran_cv_", cid, "_k10.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/ml/result/xing/marg/xing_fix_ran_cv_", cid, "_k10.rds")

    saveRDS(list(fix = t(scottFix), ran = t(scottRan), time = stime), fname1)
    saveRDS(list(fix = t(xingFix), ran = t(xingRan), time = xtime), fname2)
} else if (mtd == 6) {
    library(parallelMCMCcombine)

    cid <- id

    subJtCov <- rep(list(array(0.0, dim = c(2, 1000, 10))), 4)

    tmp <- numeric(10)

    covMap <- c(7:12, 14:18, 21:24, 28:30, 35:36, 42)
    for (kk in 1:10) {
        fname <- paste0("/Shared/ssrivastva/wasp/ml/result/comp/samp/comp_cv_", cid, "_sub_", kk, "_k10.rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        subran <- samp$samples[ , covMap]
        cov2d <- cbind(rep(2, 4), 3:6)
        for (ddd in 1:4) {
            subJtCov[[ddd]][ , , kk] <- t(subran[ , c(cov2d[ddd, 1], cov2d[ddd, 2])])
        }
        tmp[kk] <- samp$time[3]
    }

    scottCov <- list()
    scottTime <- numeric(4)
    for (ddd in 1:4) {
        strt1 <- proc.time()
        scottCov[[ddd]] <- t(consensusMCcov(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        scottTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    xingCov <- list()
    xingTime <- numeric(4)
    for (ddd in 1:4) {
        strt1 <- proc.time()
        xingCov[[ddd]] <- t(semiparamDPE(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        xingTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    fname1 <- paste0("/Shared/ssrivastva/wasp/ml/result/cons/joint/cons_cov_cv_", cid, "_k10.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/ml/result/xing/joint/xing_cov_cv_", cid, "_k10.rds")

    saveRDS(list(cov = scottCov, time = scottTime), fname1)
    saveRDS(list(cov = xingCov, time = xingTime), fname2)
} else (mtd == 7) {
    cvtrain <- readRDS("../data/ml_train.rds")
    train <- cvtrain[[id]]
    xmat <- as.matrix(train$x)
    zmat <- as.matrix(train$z)
    yvec <- as.numeric(train$y)
    group <- as.integer(train$group)

    library(inline)
    library(Rcpp)
    library(rstan)

    gg <- ordered(as.character(group), levels = sort(unique(group)))
    group <- as.integer(gg)

    simList = list(
        nobs = length(yvec),
        nfixef = ncol(xmat),
        nranef = ncol(zmat),
        ngroup = length(unique(group)),
        xmat = xmat,
        zmat = zmat,
        group = group,
        yvec = yvec)

    seeds <- (1:5000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))
    strt1 <- proc.time()
    mdl <- stan_model("full_lme.stan")
    res <- vb(mdl, data = simList, output_samples = 2000, seed = seeds[id])
    end1 <- proc.time()

    fname <- paste0("/Shared/ssrivastva/wasp/ml/result/bbvb/bbvb_res_", id, ".rds")
    saveRDS(list(res = res, time = end1 - strt1), fname)
}
