rm(list=ls())

setwd("/Shared/ssrivastva/wasp/lme/result/")

library(matrixStats)

waspCovSamp <- list()

npart <- 10

for (cc in 1:10) {
    waspCovSamp[[cc]] <- list()
    for (pp in 1:2) {
        waspCovSamp[[cc]][[pp]] <- list()
        for (kk in 1:npart) {
            cat("loaded: ", paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", pp, "_k_", kk, "_nsub10.rds"), "\n")
            dat <- readRDS(paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", pp, "_k_", kk, "_nsub10.rds"))
            cnames <- colnames(dat$samples)
            if (pp == 2) {
                waspCovSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            } else {
                waspCovSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(5:7, 9:10, 13)]
            }
        }
    }
}

for (cc in 1:10) {
    for (pp in 1:2) {
        for (kk in 1:npart) {
            for (dd in 1:3) {
                if (pp == 1) {
                    cov2d <- cbind(c(2, 2, 3), c(3, 5, 5))
                    mat <- waspCovSamp[[cc]][[pp]][[kk]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                } else {
                    cov2d <- cbind(c(2, 2, 3), c(3, 8, 8))
                    mat <- waspCovSamp[[cc]][[pp]][[kk]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                }
                cat("wrote: ", paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/joint/cov_cv_", cc, "_p_", pp, "_nsub_", kk, "_d_", dd, "_k10.csv"), "\n")
                write.table(mat, file = paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/joint/cov_cv_", cc, "_p_", pp, "_nsub_", kk, "_d_", dd, "_k10.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
            }
        }
    }
}

npart <- 20

for (cc in 1:10) {
    waspCovSamp[[cc]] <- list()
    for (pp in 1:2) {
        waspCovSamp[[cc]][[pp]] <- list()
        for (kk in 1:npart) {
            cat("loaded: ", paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", pp, "_k_", kk, "_nsub20.rds"), "\n")
            dat <- readRDS(paste0("wasp/samp/wasp_mixed_cv_", cc, "_p_", pp, "_k_", kk, "_nsub20.rds"))
            cnames <- colnames(dat$samples)
            if (pp == 2) {
                waspCovSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(81:86, 88:92, 95:98, 102:104, 109:110, 116)]
            } else {
                waspCovSamp[[cc]][[pp]][[kk]] <- dat$samples[ , c(5:7, 9:10, 13)]
            }
        }
    }
}

for (cc in 1:10) {
    for (pp in 1:2) {
        for (kk in 1:npart) {
            for (dd in 1:3) {
                if (pp == 1) {
                    cov2d <- cbind(c(2, 2, 3), c(3, 5, 5))
                    mat <- waspCovSamp[[cc]][[pp]][[kk]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                } else {
                    cov2d <- cbind(c(2, 2, 3), c(3, 8, 8))
                    mat <- waspCovSamp[[cc]][[pp]][[kk]][ , c(cov2d[dd, 1], cov2d[dd, 2])]
                }
                cat("wrote: ", paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/joint/cov_cv_", cc, "_p_", pp, "_nsub_", kk, "_d_", dd, "_k20.csv"), "\n")
                write.table(mat, file = paste0("/Shared/ssrivastva/wasp/lme/result/wasp/samp/joint/cov_cv_", cc, "_p_", pp, "_nsub_", kk, "_d_", dd, "_k20.csv"), sep = ",", row.names = FALSE, col.names = FALSE)
            }
        }
    }
}
