cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 1) {
    source("mcmc_sampler.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    cid <- cvs[id]
    did <- ndims[id]

    cvtrain <- readRDS("/Shared/ssrivastva/wasp/lme/data/mixed.rds")
    train <- cvtrain[[cid]][[did]]

    res <- sampleFromMixMdl(train$y, train$x, train$z, train$group, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/mcmc/mcmc_lme_cv_", cid, "_p_", did, ".rds")
    saveRDS(res, fname)
} else if (mtd == 2) {
    source("comp_sampler.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    tmp <- cbind(cvs, ndims)
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 10), ], rep(1:10, times = 20))

    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cid,"_k10", ".rds"))
    train <- cvtrain[[did]][[sid]]
    rm(cvtrain)

    res <- sampleFromCompMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub10.rds")
    saveRDS(res, fname)
} else if (mtd == 3) {
    source("comp_sampler.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    tmp <- cbind(cvs, ndims)
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 20), ], rep(1:20, times = 20))
    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cid,"_k20", ".rds"))
    train <- cvtrain[[did]][[sid]]
    rm(cvtrain)

    res <- sampleFromCompMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub20.rds")
    saveRDS(res, fname)
} else if (mtd == 4) {
    source("wasp_sampler.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    tmp <- cbind(cvs, ndims)
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 10), ], rep(1:10, times = 20))

    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cid,"_k10", ".rds"))
    train <- cvtrain[[did]][[sid]]
    rm(cvtrain)

    res <- sampleFromWaspMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/wasp_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub10.rds")
    saveRDS(res, fname)
} else if (mtd == 5) {
    source("wasp_sampler.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    tmp <- cbind(cvs, ndims)
    wids <- cbind(tmp[rep(1:nrow(tmp), each = 20), ], rep(1:20, times = 20))
    cid <- wids[id, 1]
    did <- wids[id, 2]
    sid <- wids[id, 3]

    cvtrain <- readRDS(paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cid,"_k20", ".rds"))
    train <- cvtrain[[did]][[sid]]
    rm(cvtrain)

    res <- sampleFromWaspMixMdl(train$y, train$x, train$z, train$group, train$nrep, 10000, 5000, 5, id)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/wasp_mixed_cv_", cid, "_p_", did, "_k_", sid, "_nsub20.rds")
    saveRDS(res, fname)
} else if (mtd == 6) {
    library(parallelMCMCcombine)
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- ndims[id]

    subfix <- array(0.0, dim = c(c(4, 80)[did], 1000, 10))
    subran <- array(0.0, dim = c(c(6, 21)[did], 1000, 10))
    tmp <- numeric(10)
    for (kk in 1:10) {
        fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", kk, "_nsub10.rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        if (did == 1) {
            subfix[ , , kk] <- t(samp$samples[ , 1:4])
            subran[ , , kk] <- t(samp$samples[ , c(5:7, 9:10, 13)])
        } else {
            subfix[ , , kk] <- t(samp$samples[ , 1:80])
            subran[ , , kk] <- t(samp$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)])
        }
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

    fname1 <- paste0("/Shared/ssrivastva/wasp/lme/result/cons/marg/cons_fix_ran_cv_", cid, "_p_", did, "_k10.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/lme/result/xing/marg/xing_fix_ran_cv_", cid, "_p_", did, "_k10.rds")

    saveRDS(list(fix = t(scottFix), ran = t(scottRan), time = stime), fname1)
    saveRDS(list(fix = t(xingFix), ran = t(xingRan), time = xtime), fname2)
} else if (mtd == 7) {
    library(parallelMCMCcombine)
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- ndims[id]

    subfix <- array(0.0, dim = c(c(4, 80)[did], 1000, 20))
    subran <- array(0.0, dim = c(c(6, 21)[did], 1000, 20))
    tmp <- numeric(20)
    for (kk in 1:20) {
        fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", kk, "_nsub20.rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        if (did == 1) {
            subfix[ , , kk] <- t(samp$samples[ , 1:4])
            subran[ , , kk] <- t(samp$samples[ , c(5:7, 9:10, 13)])
        } else {
            subfix[ , , kk] <- t(samp$samples[ , 1:80])
            subran[ , , kk] <- t(samp$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)])
        }
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

    fname1 <- paste0("/Shared/ssrivastva/wasp/lme/result/cons/marg/cons_fix_ran_cv_", cid, "_p_", did, "_k20.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/lme/result/xing/marg/xing_fix_ran_cv_", cid, "_p_", did, "_k20.rds")

    saveRDS(list(fix = t(scottFix), ran = t(scottRan), time = stime), fname1)
    saveRDS(list(fix = t(xingFix), ran = t(xingRan), time = xtime), fname2)
} else if (mtd == 8) {
    library(parallelMCMCcombine)
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- ndims[id]

    subJtCov <- rep(list(array(0.0, dim = c(2, 1000, 10))), 3)

    tmp <- numeric(10)
    for (kk in 1:10) {
        fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", kk, "_nsub10.rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        if (did == 1) {
            cov2d <- cbind(c(2, 2, 3), c(3, 5, 5))
            subcov <- samp$samples[ , c(5:7, 9:10, 13)]
            for (ddd in 1:3) {
                subJtCov[[ddd]][ , , kk] <- t(subcov[ , c(cov2d[ddd, 1], cov2d[ddd, 2])])
            }
        } else {
            cov2d <- cbind(c(2, 2, 3), c(3, 8, 8))
            subcov <- samp$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            for (ddd in 1:3) {
                subJtCov[[ddd]][ , , kk] <- t(subcov[ , c(cov2d[ddd, 1], cov2d[ddd, 2])])
            }
        }
        tmp[kk] <- samp$time[3]
    }

    scottCov <- list()
    scottTime <- numeric(3)
    for (ddd in 1:3) {
        strt1 <- proc.time()
        scottCov[[ddd]] <- t(consensusMCcov(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        scottTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    xingCov <- list()
    xingTime <- numeric(3)
    for (ddd in 1:3) {
        strt1 <- proc.time()
        xingCov[[ddd]] <- t(semiparamDPE(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        xingTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    fname1 <- paste0("/Shared/ssrivastva/wasp/lme/result/cons/joint/cons_cov_cv_", cid, "_p_", did, "_k10.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/lme/result/xing/joint/xing_cov_cv_", cid, "_p_", did, "_k10.rds")

    saveRDS(list(cov = scottCov, time = scottTime), fname1)
    saveRDS(list(cov = xingCov, time = xingTime), fname2)
} else if (mtd == 9) {
    library(parallelMCMCcombine)
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)
    cid <- cvs[id]
    did <- ndims[id]

    subJtCov <- rep(list(array(0.0, dim = c(2, 1000, 20))), 3)

    tmp <- numeric(20)
    for (kk in 1:20) {
        fname <- paste0("/Shared/ssrivastva/wasp/lme/result/comp/comp_mixed_cv_", cid, "_p_", did, "_k_", kk, "_nsub20.rds")
        samp <- readRDS(fname)
        cnames <- colnames(samp$samples)
        if (did == 1) {
            cov2d <- cbind(c(2, 2, 3), c(3, 5, 5))
            subcov <- samp$samples[ , c(5:7, 9:10, 13)]
            for (ddd in 1:3) {
                subJtCov[[ddd]][ , , kk] <- t(subcov[ , c(cov2d[ddd, 1], cov2d[ddd, 2])])
            }
        } else {
            cov2d <- cbind(c(2, 2, 3), c(3, 8, 8))
            subcov <- samp$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            for (ddd in 1:3) {
                subJtCov[[ddd]][ , , kk] <- t(subcov[ , c(cov2d[ddd, 1], cov2d[ddd, 2])])
            }
        }
        tmp[kk] <- samp$time[3]
    }

    scottCov <- list()
    scottTime <- numeric(3)
    for (ddd in 1:3) {
        strt1 <- proc.time()
        scottCov[[ddd]] <- t(consensusMCcov(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        scottTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    xingCov <- list()
    xingTime <- numeric(3)
    for (ddd in 1:3) {
        strt1 <- proc.time()
        xingCov[[ddd]] <- t(semiparamDPE(subchain = subJtCov[[ddd]]))
        end1 <- proc.time()
        xingTime[ddd] <- mean(tmp) + end1[3] - strt1[3]
    }

    fname1 <- paste0("/Shared/ssrivastva/wasp/lme/result/cons/joint/cons_cov_cv_", cid, "_p_", did, "_k20.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/lme/result/xing/joint/xing_cov_cv_", cid, "_p_", did, "_k20.rds")

    saveRDS(list(cov = scottCov, time = scottTime), fname1)
    saveRDS(list(cov = xingCov, time = xingTime), fname2)
} else if (mtd == 10) {
    source("variational_bayes.R")
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    cid <- cvs[id]
    did <- ndims[id]

    cvtrain <- readRDS("/Shared/ssrivastva/wasp/lme/data/mixed.rds")
    train <- cvtrain[[cid]][[did]]

    res <- fitLinearMixefEffectsVB(train$y, train$x, train$z, train$group, 1000)
    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/vb/vb_lme_cv_", cid, "_p_", did, ".rds")
    saveRDS(res, fname)
} else (mtd == 11) {
    cvs <- rep(1:10, each = 2)
    ndims <- rep(1:2, times = 10)

    cid <- cvs[id]
    did <- ndims[id]

    cvtrain <- readRDS("/Shared/ssrivastva/wasp/lme/data/mixed.rds")
    train <- cvtrain[[cid]][[did]]

    library(inline)
    library(Rcpp)
    library(rstan)

    yvec <- train$y; xmat <- train$x; zmat <- train$z; group <- train$group
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

    fname <- paste0("/Shared/ssrivastva/wasp/lme/result/bbvb/bbvb_lme_cv_", cid, "_p_", did, ".rds")
    saveRDS(list(res = res, time = end1 - strt1), fname)
}
