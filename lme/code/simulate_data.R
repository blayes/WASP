rm(list=ls())

set.seed(12345)

setwd("~/wasp/lme/code")
library(matrixStats)
library(Matrix)

## See mbest package at http://ptrckprry.com/code/ Perry (2017) in JRSS-B.
genData <- function (ngroup, nobs, nfixef, nranef) {
    ## fixed effects coefficients
    fixef <- rep(c(-2, 2), length = nfixef)
    if (nranef == 3) {
        ranefCorr <- matrix(c(1, -0.4, 0.3,
                              -0.4, 1, 0.001,
                              0.3, 0.001, 1),
                            nranef, nranef)
    } else {
        ranefCorr <- as.matrix(bdiag(rep(list(matrix(c(1, -0.4, 0.3,
                                                       -0.4, 1, 0.001,
                                                       0.3, 0.001, 1),
                                                     3, 3)), 2)))
    }
    ranefCov <- outer(sqrt(1:nranef), sqrt(1:nranef)) * ranefCorr
    ranefCovSqrt <- chol(ranefCov)

    # generate coefficients
    u <- matrix(rnorm(ngroup * nranef), ngroup, nranef)
    ranef <- u %*% ranefCovSqrt

    ## generate group
    suppressWarnings({ # ignore warning about using Walker's alias method
        group <- sample.int(ngroup, nobs, replace=TRUE)
    })

    ## generate feature  matrices with Pr(x[i,j] = +1) = P(x[i,j] = -1) = 1/2,
    x <- matrix(sample(c(-1, +1), nobs * nfixef, replace=TRUE), nobs, nfixef)
    z <- matrix(sample(c(-1, +1), nobs * nranef, replace=TRUE), nobs, nranef)

    ## compute linear predictors and generate observations
    mu <- drop(x %*% fixef) + rowSums(z * ranef[group,])
    y <- rnorm(nobs, mean=mu, sd=1)

    list(ngroup = ngroup, nobs = nobs,
         #fixef = fixef, #ranef = ranef,
         ranefCov = ranefCov,
         ranefCovSqrt = ranefCovSqrt,
         group = group, x = x, z = z, y.mean = mu, y = y)
}

ngroup <- 6000
nobs <- 1e5
nfixef <- c(4, 80)
nranef <- c(3, 6)

repData <- list()
for (cc in 1:10) {
  repData[[cc]] <- vector("list", 2)
  names(repData[[cc]]) <- paste0("p", nfixef, "q", nranef)
}

for (cc in 1:10) {
  cat("cc ", cc, "\n")
  for (pp in 1:2) {
      repData[[cc]][[pp]] <- genData(ngroup, nobs, nfixef[pp], nranef[pp])
  }
}

saveRDS(repData, "/Shared/ssrivastva/wasp/lme/data/mixed.rds")

rm(list=ls())

repData <- readRDS("/Shared/ssrivastva/wasp/lme/data/mixed.rds")

set.seed(12345)

ngroup <- 6000
nobs <- 1e5
nfixef <- c(4, 80)
nranef <- c(3, 6)

npart <- 10
partData <- list()

for (cc in 1:10) {
    partData <- vector("list", 2)
    names(partData) <- paste0("p", nfixef, "q", nranef)
    for (pp in 1:2) {
        partData[[pp]] <- vector("list", npart)
        names(partData[[pp]]) <- paste0("k", 1:npart)
        lst <- repData[[cc]][[pp]]
        grpSplit <- split(1:nrow(lst$x), lst$group)
        partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
        for (ll in 1:npart) {
            grpIdx <- which(partsIdx == ll)
            idx <- unlist(grpSplit[grpIdx])
            partData[[pp]][[ll]]$nobs <- length(idx)
            partData[[pp]][[ll]]$x <- lst$x[idx, ]
            partData[[pp]][[ll]]$y <- lst$y[idx]
            partData[[pp]][[ll]]$z <- lst$z[idx, ]
            partData[[pp]][[ll]]$group <- lst$group[idx]
            partData[[pp]][[ll]]$idx <- idx
            partData[[pp]][[ll]]$nrep <- nobs / length(idx)
        }
    }
    saveRDS(partData, paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cc, "_k10", ".rds"))
}

rm(list=ls())

repData <- readRDS("/Shared/ssrivastva/wasp/lme/data/mixed.rds")

set.seed(12345)

ngroup <- 6000
nobs <- 1e5
nfixef <- c(4, 80)
nranef <- c(3, 6)

npart <- 20
partData <- list()

for (cc in 1:10) {
    partData <- vector("list", 2)
    names(partData) <- paste0("p", nfixef, "q", nranef)
    for (pp in 1:2) {
        partData[[pp]] <- vector("list", npart)
        names(partData[[pp]]) <- paste0("k", 1:npart)
        lst <- repData[[cc]][[pp]]
        grpSplit <- split(1:nrow(lst$x), lst$group)
        partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
        for (ll in 1:npart) {
            grpIdx <- which(partsIdx == ll)
            idx <- unlist(grpSplit[grpIdx])
            partData[[pp]][[ll]]$nobs <- length(idx)
            partData[[pp]][[ll]]$x <- lst$x[idx, ]
            partData[[pp]][[ll]]$y <- lst$y[idx]
            partData[[pp]][[ll]]$z <- lst$z[idx, ]
            partData[[pp]][[ll]]$group <- lst$group[idx]
            partData[[pp]][[ll]]$idx <- idx
            partData[[pp]][[ll]]$nrep <- nobs / length(idx)
        }
    }
    saveRDS(partData, paste0("/Shared/ssrivastva/wasp/lme/data/wasp_mixed_cv_", cc, "_k20", ".rds"))
    cat("cc: ", cc, "\n")
}
