rm(list=ls())
setwd("/Shared/ssrivastva/wasp/ml/result/")

meanMap <- 1:6
covMap <- c(7:12, 14:18, 21:24, 28:30, 35:36, 42)

meanList <- list()
covList <- list()
for (cc in 1:10) {
    meanList[[cc]] <- list()
    covList[[cc]] <- list()
    for (kk in 1:10) {
        meanList[[cc]][[kk]] <- readRDS(paste0("wasp/samp/wasp_cv_", cc, "_sub_", kk, "_k10.rds"))$samples[ , meanMap]
        covList[[cc]][[kk]] <- readRDS(paste0("wasp/samp/wasp_cv_", cc, "_sub_", kk, "_k10.rds"))$samples[ , covMap]
    }
}

cov2d <- cbind(rep(2, 4), 3:6)

for (cc in 1:10) {
    for (kk in 1:10) {
        for (ddd in 1:4) {
            write.table(covList[[cc]][[kk]][ , c(cov2d[ddd, 1], cov2d[ddd, 2])],
                        file = paste0("wasp/samp/joint/cov_cv_", cc, "_nsub_", kk, "_d_", ddd,"_k10.csv"),
                        sep = ",", row.names = FALSE, col.names = FALSE)
        }
    }
}
