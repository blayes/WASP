rm(list = ls())
setwd("~/wasp/ml/code/")
fdata <- readRDS("~/wasp/ml/data/dataSel.rds")

set.seed(12345)

ncv <- 10
cvgrps <- sample(1:ncv, length(unique(fdata$group)), replace = TRUE)
grpSplit <- split(1:nrow(fdata$x), fdata$group)

testIdx <- list()
trainIdx <- list()
train <- list()
test <- list()
for (ii in 1:10) {
    idx <- (cvgrps == ii)
    testIdx[[ii]] <- unlist(grpSplit[idx])
    test[[ii]] <- list(x = fdata$x[unlist(grpSplit[idx]), ],
                       z = fdata$z[unlist(grpSplit[idx]), ],
                       y = fdata$y[unlist(grpSplit[idx])],
                       group = fdata$group[unlist(grpSplit[idx])]
                       )
    trainIdx[[ii]] <- setdiff(1:nrow(fdata$x), unlist(grpSplit[idx]))
    train[[ii]] <- list(x = fdata$x[-unlist(grpSplit[idx]), ],
                        z = fdata$z[-unlist(grpSplit[idx]), ],
                        y = fdata$y[-unlist(grpSplit[idx])],
                        group = fdata$group[-unlist(grpSplit[idx])]
                        )
}

saveRDS(train, "../data/ml_train.rds")
saveRDS(test, "../data/ml_test.rds")

saveRDS(testIdx, "../data/test_idx.rds")
saveRDS(trainIdx, "../data/train_idx.rds")

rm(list = ls())

train <- readRDS("../data/ml_train.rds")

ncv <- 10
npart <- 10
parts <- list()
for (jj in 1:ncv) {
    parts[[jj]] <- vector("list", length = npart)
    names(parts[[jj]]) <- paste0("part", 1:npart)
}
names(parts) <- paste0("cv", 1:ncv)

set.seed(12345)
npart <- 10
parts <- list()
for (cc in 1:10) {
    parts[[cc]] <- vector("list", npart)
    names(parts[[cc]]) <- paste0("k", 1:npart)
    lst <- train[[cc]]
    grpSplit <- split(1:nrow(lst$x), lst$group)
    partsIdx <- sample(1:npart, length(grpSplit), replace = TRUE)
    for (ll in 1:npart) {
        grpIdx <- which(partsIdx == ll)
        idx <- unlist(grpSplit[grpIdx])
        parts[[cc]][[ll]]$nobs <- length(idx)
        parts[[cc]][[ll]]$x <- lst$x[idx, ]
        parts[[cc]][[ll]]$y <- lst$y[idx]
        parts[[cc]][[ll]]$z <- lst$z[idx, ]
        parts[[cc]][[ll]]$group <- lst$group[idx]
        parts[[cc]][[ll]]$idx <- idx
        parts[[cc]][[ll]]$nrep <- nrow(train[[cc]]$x) / length(idx)
    }
}

saveRDS(parts, "../data/wasp_ml_train.rds")
