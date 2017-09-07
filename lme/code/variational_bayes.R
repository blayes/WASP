fitLinearMixefEffectsVB <- function (yvec, xmat, zmat, group, niter) {
    library(Matrix)
    library(MCMCpack)

    nfixef <- ncol(xmat)
    nranef <- ncol(zmat)
    grpLbl <- sort(unique(group))
    ngroup <- length(grpLbl)
    ndim <- nrow(xmat)

    elbo <- numeric(niter)

    ranefList <- list()
    grpIdx <- list()
    for (ii in 1:ngroup) {
        grpIdx[[ii]] <- which(group == grpLbl[ii])
        ranefList[[ii]] <- zmat[grpIdx[[ii]], , drop = FALSE]
    }
    ranefMat <- bdiag(ranefList)

    fixefMat <- xmat[unlist(grpIdx), ]

    designMat <- cBind(fixefMat, ranefMat)
    designTransDesign <- crossprod(designMat, designMat)
    ordY <- yvec[unlist(grpIdx)]
    ordGrp <- group[unlist(grpIdx)]
    designTransY <- crossprod(designMat, ordY)

    muErrInv <- 1; muErrAInv <- 1; muVarErrAinv <- 1; errAScl <- 1; nu = 2; muRanCovInv <- solve(rWishart(1, 2 * ncol(zmat), diag(ncol(zmat)))[ , , 1]); varBeta <- 100; ranCovAsScl <- rep(1, ncol(zmat));

    prev <- list(coefMu = runif(nfixef), coefCov = diag(nfixef), ranCov = diag(nranef)); conv <- 1e5

    startTime <- proc.time()
    for (its in 0:niter) {

        smat <- diag(0, ncol(xmat))
        svec <- numeric(ncol(xmat))

        gmat <- vector("list", ngroup)
        hmat <- vector("list", ngroup)
        for(ii in 1:ngroup) {
            gmat[[ii]] <- crossprod(xmat[grpIdx[[ii]], ], zmat[grpIdx[[ii]], ]) * muErrInv
            hmat[[ii]] <- solve(crossprod(zmat[grpIdx[[ii]], ], zmat[grpIdx[[ii]], ]) * muErrInv + muRanCovInv)
            tmp <- gmat[[ii]] %*% hmat[[ii]]
            smat <- smat + tcrossprod(tmp, gmat[[ii]])
            svec <- svec + drop(tmp %*% crossprod(zmat[grpIdx[[ii]], ], yvec[grpIdx[[ii]]]))
        }

        fixCov <- solve(crossprod(fixefMat, fixefMat) * muErrInv + diag(1, nfixef) / varBeta - smat)
        fixMu <- drop(muErrInv * fixCov %*% (crossprod(fixefMat, ordY) - svec))

        ranMu <- vector("list", ngroup)
        ranCov <- vector("list", ngroup)
        tmp1 <- numeric(ngroup)
        tmp2 <- numeric(ngroup)
        for(ii in 1:ngroup) {
            ranCov[[ii]] <- hmat[[ii]] + tcrossprod(hmat[[ii]], gmat[[ii]]) %*% fixCov %*% gmat[[ii]] %*% hmat[[ii]]
            ranMu[[ii]] <- drop(hmat[[ii]] %*% (muErrInv * crossprod(zmat[grpIdx[[ii]], ], yvec[grpIdx[[ii]]]) - crossprod(gmat[[ii]], fixMu)))
            tmp1[ii] <- sum(crossprod(zmat[grpIdx[[ii]], ], zmat[grpIdx[[ii]], ]) * ranCov[[ii]])
            tmp2[ii] <- sum(tcrossprod(gmat[[ii]] %*% hmat[[ii]], gmat[[ii]]) * fixCov)
        }
        ranMuVec <- do.call(c, ranMu)
        resids <- ordY - fixefMat %*% fixMu - ranefMat %*% ranMuVec
        scaleErr <- muErrAInv + 0.5 * (sum(resids^2) + sum(crossprod(fixefMat, fixefMat) * fixCov) + sum(tmp1) - 2 * sum(tmp2) / muErrInv)
        shapeErr <- 0.5 * (ndim + 1)
        muErrInv <- shapeErr / scaleErr

        ## update post. for parameter in the px-ed form of half-cauchy prior for err
        shapeErrA <- 1
        scaleErrA <- muErrInv + errAScl
        muErrAInv <- shapeErrA / scaleErrA

        ## update post. for parameter in the px-ed form of half-cauchy
        ## prior for the random effects covariance matrix
        scaleRanCovAs <- nu * diag(muRanCovInv) + ranCovAsScl
        shapeRanCovAs <- 0.5 * (nu + nranef)
        muRanCovAsInv <- pmax(shapeRanCovAs / scaleRanCovAs, 1e-5)

        ranefCoefMuMat <- matrix(ranMuVec, nrow = nranef, ncol = ngroup)
        ranefCoefCovMat <- diag(0, nranef)
        for (jj in 1:nranef) {
            ranefCoefCovMat <- ranefCoefCovMat + ranCov[[jj]]
        }
        scaleRanCovMat <- tcrossprod(ranefCoefMuMat, ranefCoefMuMat) + ranefCoefCovMat + 2 * nu * diag(muRanCovAsInv)
        rateRanCovMat <- solve(scaleRanCovMat)
        muRanCovInv <- (nu + ngroup + nranef - 1) * rateRanCovMat

        if ((its > 10) && (conv < 1e-10))
            break

        if (its %% 10 == 0) {
            cat("iteration: ", its, "\n")

            diff1 <- matrix(prev$coefMu[1:nfixef] - fixMu)
            diff2 <- prev$coefCov - fixCov
            diff3 <- matrix(prev$ranCov - as.matrix(scaleRanCovMat))
            conv <- norm(diff1, "O") + norm(diff2, "O") + norm(diff3, "O")

            prev$coefMu <- fixMu; prev$coefCov <- fixCov; prev$ranCov <- scaleRanCovMat
        }

    }
    endTime <- proc.time()

    list(
        coefs = list(
            cov = fixCov,
            mu = fixMu
            )
        ,
        err = list(
            aa = scaleErr,
            bb = shapeErr
            )
        ,
        cov = list(
            scale = scaleRanCovMat,
            df = (nu + ngroup + nranef - 1)
            )
        ,
        niter = its
        ,
        time = endTime - startTime
        )
}
