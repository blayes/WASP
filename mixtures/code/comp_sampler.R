mvnCompMix <- function (dataMat, ncomp = 2, nrep, niter = 10000, nburn = 5000, nthin = 5) {
    library(matrixStats)
    library(mvtnorm)
    library(MCMCpack)

    nobs <- nrow(dataMat)
    ndim <- ncol(dataMat)
    alpha <- rep(1 / ncomp, ncomp)

    probs <- as.numeric(rdirichlet(1, alpha))
    zMat <- rmultinom(nobs, 1, probs)
    densMat <- matrix(0.0, nrow = nobs, ncol = ncomp)
    muMat <- matrix(0.0, nrow = ncomp, ncol = ndim)
    sigArr <- aperm(array(diag(1.0, ndim), c(ndim, ndim, ncomp)), perm = c(3, 1, 2))

    kap0 <- 0.01; m0 <- rep(0.0, ndim) # mu hyper-pars
    nu0 <- 2; s0 <- 2 * nu0 * diag(1.0, ndim) # sigma hyper-pars

    probsSamp <- matrix(0.0, nrow = (niter - nburn) / nthin, ncol = ncomp)
    muMatSamp <- array(0.0, dim = c(ncomp, ndim, (niter - nburn) / nthin))
    sigMatSamp <- array(0.0, dim = c(ncomp, ndim, ndim, (niter - nburn) / nthin))

    cts <- 0
    startTime <- proc.time()
    for (ii in 0:niter) {
        idxList <- lapply(split(zMat, row(zMat)),
                          function (x) {
                              which(x == 1)
                          })
        ns <- sapply(idxList, length)
        # sample probs
        probs <- rdirichlet(1, ns + (alpha - 1 + nrep) / nrep)

        if (ii %% 100 == 0) {cat("gibbs: ", ii, "\n"); cat("ns: ", ns, "\n")}

        for (jj in seq_along(ns)) {
            datMean <- colMeans(dataMat[idxList[[jj]], , drop = FALSE])
            ## mean
            muCov <- sigArr[jj, , ] / (kap0 / nrep + ns[jj])
            muMean <- (kap0 * m0 / nrep + ns[jj] * datMean) / (kap0 / nrep + ns[jj])
            muMat[jj, ] <- rmvnorm(1, mean = muMean, sigma = muCov, method = "chol")
            ## cov
            mat1 <- (kap0 / nrep * ns[jj] /  (kap0 / nrep + ns[jj])) * tcrossprod(datMean - m0, datMean - m0)
            centMat <- dataMat[idxList[[jj]], , drop = FALSE] - matrix(datMean, nrow = length(idxList[[jj]]), ncol = ndim, byrow = TRUE)
            mat2 <- crossprod(centMat, centMat)
            covSclMat <- mat1 + mat2 + s0 / nrep
            covDf <- ns[jj] + (nu0 + ndim + 1) / nrep - (ndim + 1)
            sigArr[jj, , ]<- riwish(covDf, covSclMat)
            densMat[ , jj] <- dmvnorm(dataMat, muMat[jj, ], sigArr[jj, , ], log = TRUE)
        }

        lprobs <- densMat + log(matrix(probs, nrow = nobs, ncol = ncol(densMat), byrow = TRUE))
        eprobs <- exp(lprobs - rowMaxs(lprobs)) / rowSums(exp(lprobs - rowMaxs(lprobs)))

        for (kk in seq_len(nobs)) {
            ppp <- eprobs[kk, ]
            zMat[ , kk] <- rmultinom(1, 1, ppp)
        }

        if ((ii > nburn) && (ii %% nthin == 0)) {
            cts <- cts + 1
            probsSamp[cts, ] <- probs
            muMatSamp[ , , cts] <- muMat
            sigMatSamp[ , , , cts] <- sigArr
        }
    }
    endTime <- proc.time()


    list(
        'mu' = muMatSamp,
        'cov' = sigMatSamp,
        'prob' = probsSamp,
        'time' = endTime[3] - startTime[3]
    )
}
