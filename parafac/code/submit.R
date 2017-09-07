cmdArgs <- commandArgs(trailingOnly = TRUE)

mtd <- as.numeric(cmdArgs[1])
id <- as.numeric(cmdArgs[2])

if (mtd == 2) {
    library(parallelMCMCcombine)

    cid <- id
    subProb <- array(0.0, dim = c(20, 1000, 5))

    for (dd in 1:20) {
        for (kk in 1:5) {
            dat <- read.table(paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub5/res_cv_", cid, "_sub_", kk, "_dim_", dd, "_k5.csv"), sep = ",", header = FALSE)
            subProb[dd, , kk] <- dat[ , 1]
        }
    }

    stime <- rep(NA, 20)
    scotRes <- list()
    for (dd in 1:20) {
        strt1 <- proc.time()
        scotRes[[dd]] <- as.numeric(t(consensusMCindep(subchain = subProb[dd, , , drop = FALSE])))
        end1 <- proc.time()
        stime[dd] <- end1[3] - strt1[3]
    }

    xtime <- rep(NA, 20)
    xingRes <- list()
    for (dd in 1:20) {
        strt1 <- proc.time()
        xingRes[[dd]] <- as.numeric(t(semiparamDPE(subchain = subProb[dd, , , drop = FALSE])))
        end1 <- proc.time()
        xtime[dd] <- end1[3] - strt1[3]
    }

    fname1 <- paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub5/marg/cons_cv_", cid, "_k5.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub5/marg/xing_cv_", cid, "_k5.rds")

    saveRDS(list(marg = scotRes, time = stime), fname1)
    saveRDS(list(marg = xingRes, time = xtime), fname2)
} else if (mtd == 3) {
    library(parallelMCMCcombine)

    cid <- id
    subProb <- array(0.0, dim = c(20, 1000, 10))

    for (dd in 1:20) {
        for (kk in 1:10) {
            dat <- read.table(paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub10/res_cv_", cid, "_sub_", kk, "_dim_", dd, "_k10.csv"), sep = ",", header = FALSE)
            subProb[dd, , kk] <- dat[ , 1]
        }
    }

    stime <- rep(NA, 20)
    scotRes <- list()
    for (dd in 1:20) {
        strt1 <- proc.time()
        scotRes[[dd]] <- as.numeric(t(consensusMCindep(subchain = subProb[dd, , , drop = FALSE])))
        end1 <- proc.time()
        stime[dd] <- end1[3] - strt1[3]
    }

    xtime <- rep(NA, 20)
    xingRes <- list()
    for (dd in 1:20) {
        strt1 <- proc.time()
        xingRes[[dd]] <- as.numeric(t(semiparamDPE(subchain = subProb[dd, , , drop = FALSE])))
        end1 <- proc.time()
        xtime[dd] <- end1[3] - strt1[3]
    }

    fname1 <- paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub10/marg/cons_cv_", cid, "_k10.rds")
    fname2 <- paste0("/Shared/ssrivastva/wasp/parafac/result/comp/sub10/marg/xing_cv_", cid, "_k10.rds")

    saveRDS(list(marg = scotRes, time = stime), fname1)
    saveRDS(list(marg = xingRes, time = xtime), fname2)
}
