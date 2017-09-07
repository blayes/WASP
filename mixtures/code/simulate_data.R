rm(list=ls())
setwd("/Users/ssrivastva/wasp/mixtures/code/")

genData <- function (muList = list('1' = c(1, 2), '2' = c(7, 8)),
                     sigMat = matrix(c(1, 0.5, 0.5, 2), 2, 2),
                     probs = c(0.3, 0.7),
                     nobs = 1000) {
    library(mvtnorm)
    zs <- rmultinom(nobs, 1, probs)

    idxList <- lapply(split(zs, row(zs)),
                      function (x) {
                          which(x == 1)
                      })



    dataMat <- matrix(0.0, nrow = nobs, ncol = 2)
    clusts <- numeric(nobs)
    for (ii in seq_along(muList)) {
        idx <- idxList[[ii]]
        clusts[idx] <- ii
        dataMat[idx, ] <- rmvnorm(length(idx), muList[[ii]], sigMat)
    }
    rownames(dataMat) <- paste("data", seq_len(nrow(dataMat)), clusts, sep = "_")

    list(data = dataMat,
         cluster = clusts,
         zs = zs
         )
}

#### full data

set.seed(12345)
nrep <- 10
reps <- vector("list", length = nrep)
nobs <- 1e5
for (r in seq_len(nrep)) {
    cat(r, "rep\n")
    reps[[r]] <- genData(nobs = nobs)
}
saveRDS(reps, "../data/mix_100k.rds")

#### k = 10

rm(list = ls())
reps <- readRDS("../data/mix_100k.rds")

npart <- 10
nrep <- 10
nclust <- 3

parts <- vector("list", length = nrep)
partsIdx <- vector("list", length = nrep)

for (ii in 1:nrep) {
    parts[[ii]] <- vector("list", length = npart)
}

set.seed(12345)
for (r in seq_len(nrep)) {
    kmns <- kmeans(reps[[r]]$data, nclust)
    partsIdx[[r]] <- numeric(nrow(reps[[r]]$data))

    for (cc in 1:nclust) {
        ccIdx <- which(kmns$cluster == cc)
        partsIdx[[r]][ccIdx] <- sample(1:npart, length(ccIdx), replace = TRUE)
    }

    for (ii in 1:npart) {
        parts[[r]][[ii]] <- reps[[r]]$data[partsIdx[[r]] == ii, ]
    }
}

saveRDS(parts, "../data/wasp_mix_100k_k10.rds")

#### k = 5

rm(list = ls())
reps <- readRDS("../data/mix_100k.rds")

npart <- 5
nrep <- 10
nclust <- 3

parts <- vector("list", length = nrep)
partsIdx <- vector("list", length = nrep)

for (ii in 1:nrep) {
    parts[[ii]] <- vector("list", length = npart)
}

set.seed(12345)
for (r in seq_len(nrep)) {
    kmns <- kmeans(reps[[r]]$data, nclust)
    partsIdx[[r]] <- numeric(nrow(reps[[r]]$data))

    for (cc in 1:nclust) {
        ccIdx <- which(kmns$cluster == cc)
        partsIdx[[r]][ccIdx] <- sample(1:npart, length(ccIdx), replace = TRUE)
    }

    for (ii in 1:npart) {
        parts[[r]][[ii]] <- reps[[r]]$data[partsIdx[[r]] == ii, ]
    }
}

saveRDS(parts, "../data/wasp_mix_100k_k5.rds")
