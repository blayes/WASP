sampleFromCompMixMdl <- function (yvec, xmat, zmat, group, nrep, niter, nburn, nthin, id) {
    library(inline)
    library(Rcpp)
    library(rstan)

    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)

    ranefList <- list()
    grpIdx <- list()
    for (ii in 1:ngroup) {
      grpIdx[[ii]] <- which(group == grpLbl[ii])
      ranefList[[ii]] <- zmat[grpIdx[[ii]], , drop = FALSE]
    }
    ranefMat <- do.call(rbind, ranefList)
    fixefMat <- xmat[unlist(grpIdx), ]

    pos2 <- cumsum(sapply(grpIdx, length))
    pos1 <- c(1, pos2[-ngroup] + 1)

    ordY <- yvec[unlist(grpIdx)]
    ordGrp <- group[unlist(grpIdx)]

    idx <- seq_along(unlist(grpIdx))
    simList = list(nobs = length(yvec[idx]),
      nfixef = ncol(xmat),
      nranef = ncol(zmat),
      ngroup = length(unique(ordGrp[idx])),
      nsub = nrep,
      xmat = fixefMat[idx, ],
      zmat = ranefMat[idx, ],
      group = ordGrp[idx],
      yvec = ordY[idx],
      pos1 = pos1,
      pos2 = pos2)

    seeds <- (1:2000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))

    startTime <- proc.time()
    mdl <- stan(file = "comp_lme.stan", data = simList, iter = niter, warmup = nburn, chains = 1, thin = nthin,
                seed = seeds[id],
                init = list(list(betas = rep(0, ncol(xmat)),
                                 corrRanef = diag(ncol(zmat)),
                                 sclRanef = rep(2, ncol(zmat))
                                 )))
    endTime <- proc.time()

    lst <- mdl@sim$samples[[1]]
    bs <- grep("fixef|covRanef", names(lst))
    sampdf <- do.call(cbind, lst[bs])

    list(samples = sampdf[(nrow(sampdf) - (niter - nburn) / nthin + 1):nrow(sampdf), ], time = endTime - startTime)
}
