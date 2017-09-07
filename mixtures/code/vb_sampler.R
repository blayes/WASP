# Bishop's approach copied from http://www.cs.ubc.ca/~murphyk/Software/VBEMGMM/index.html
mvnVbMix <- function (dataMat, ncomp = 2, niter = 1000) {
    library(matrixStats)
    library(mvtnorm)
    library(MCMCpack)

    nobs <- nrow(dataMat)
    ndim <- ncol(dataMat)

    alpha0 <- rep(1 / ncomp, ncomp)
    beta0 <- 1; m0 <- rep(0.0, ndim) # mu hyper-pars;
    v0 <- 2; W0inv <- 2 * v0 * diag(1.0, ndim); W0 <- solve(W0inv) # sigma hyper-pars

    rho <- rdirichlet(nobs, alpha0)
    r <- rho
    Nk <- colSums(r)
    xbar <- vector("list", ncomp)
    S <- vector("list", ncomp)
    for (kk in 1:ncomp) {
        xbar[[kk]] <- colSums(r[ , kk] * dataMat) / Nk[kk]
        S[[kk]] <- crossprod((dataMat - matrix(xbar[[kk]], nrow = nobs, ncol = ndim, byrow = TRUE)) * r[ , kk], (dataMat - matrix(xbar[[kk]], nrow = nobs, ncol = ndim, byrow = TRUE))) / Nk[kk]
    }

    alpha <- alpha0 + Nk
    betas <- beta0 + Nk
    v <- v0 + Nk
    m <- vector("list", ncomp)
    W <- vector("list", ncomp)
    for (kk in 1:ncomp) {
        m[[kk]] <- (beta0 * m0 + Nk[kk] * xbar[[kk]]) / betas[kk]
        mlt <- (beta0 * Nk[kk]) / (beta0 + Nk[kk])
        W[[kk]] <- solve(W0inv + Nk[kk] * S[[kk]] + mlt * tcrossprod(xbar[[kk]] - m0))
    }

    logLambdaTilde <- numeric(ncomp)
    E <- matrix(NA, nobs, ncomp)
    trSW <- numeric(ncomp)
    xbarWxbar <- numeric(ncomp)
    mWm <- numeric(ncomp)
    trW0invW <- numeric(ncomp)
    L <- numeric(niter + 1)

    startTime <- proc.time()
    for (ii in 0:niter) {
        if (ii %% 10 == 0) cat("iter: ", ii, "\n")

        psiAlphaHat <- psigamma(sum(alpha))
        logPiTilde <- psigamma(alpha) - psiAlphaHat
        const <- ndim * log(2)
        for (kk in 1:ncomp) {
            logLambdaTilde[kk] <- sum(psigamma((v[kk] - (1:ndim) + 1) * 0.5)) + const + as.numeric(determinant(W[[kk]])$modulus)
            for (nn in 1:nobs) {
                tmp1 <- dataMat[nn, ] - m[[kk]]
                E[nn, kk] <- (ndim / betas[kk]) + v[kk] * crossprod(tmp1, W[[kk]] %*% tmp1)
            }
        }

        logRho <- matrix(logPiTilde, nobs, ncomp, byrow = TRUE) + 0.5 * matrix(logLambdaTilde, nobs, ncomp, byrow = TRUE) - 0.5 * E
        r <- exp(logRho - rowMaxs(logRho)) / rowSums(exp(logRho - rowMaxs(logRho)))

        Nk <- colSums(r)
        for (kk in 1:ncomp) {
            xbar[[kk]] <- colSums(r[ , kk] * dataMat) / Nk[kk]
            S[[kk]] <- crossprod((dataMat - matrix(xbar[[kk]], nrow = nobs, ncol = ndim, byrow = TRUE)) * r[ , kk], (dataMat - matrix(xbar[[kk]], nrow = nobs, ncol = ndim, byrow = TRUE))) / Nk[kk]
        }

        logCalpha0 <- lgamma(ncomp * alpha0) - ncomp * lgamma(alpha0)
        logB0 <- (v0 / 2) * as.numeric(determinant(W0inv)$modulus) - (v0 * ndim/2) * log(2) - (ndim*(ndim-1)/4) * log(pi) - sum(lgamma(0.5 * (v0 + 1 - (1:ndim))))
        logCalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))

        H <- 0
        for (kk in 1:ncomp) {
            logBk <- -(v[kk]/2)*log(det(W[[kk]])) - (v[kk]*ndim/2)*log(2) - (ndim*(ndim-1)/4)*log(pi) - sum(lgamma(0.5*(v[kk] + 1 - (1:ndim))))
            H <- H -logBk - 0.5*(v[kk] -ndim-1) * logLambdaTilde[kk] + 0.5*v[kk]*ndim
            trSW[kk] <- sum(diag(v[kk] * S[[kk]] %*% W[[kk]]))
            diff1 <- xbar[[kk]] - m[[kk]]
            xbarWxbar[[kk]] <- sum(diff1 * (W[[kk]] %*% diff1))
            diff1 <- m[[kk]] - m0
            mWm[kk] <- sum(diff1 * (W[[kk]] %*% diff1))
            trW0invW[kk] <- sum(diag(W0inv %*% W[[kk]]))
        }

        Lt1 <- 0.5 * sum(Nk * (logLambdaTilde - ndim / betas - trSW - v * xbarWxbar - ndim * log(2*pi)));
        Lt2 <- sum(Nk * logPiTilde)
        Lt3 <- sum(logCalpha0) + sum((alpha0 - 1) * sum(logPiTilde))
        Lt41 <- 0.5 * sum(ndim * log(beta0/(2*pi)) + logLambdaTilde - ndim * beta0 / betas - beta0 * v * mWm)
        Lt42 <- ncomp * logB0 + 0.5 * (v0-ndim-1) * sum(logLambdaTilde) - 0.5*sum(v * trW0invW)
        Lt4 <- Lt41 + Lt42
        Lt5 <- sum(sum(r * log(r)))
        Lt6 <- sum((alpha - 1) * logPiTilde) + logCalpha
        Lt7 <- 0.5 * sum(logLambdaTilde + ndim * log(betas / (2*pi))) - 0.5 * ndim * ncomp - H

        L[ii] <- Lt1 + Lt2 + Lt3 + Lt4 - Lt5 - Lt6 - Lt7

        alpha <- alpha0 + Nk
        betas <- beta0 + Nk
        v <- v0 + Nk
        m <- vector("list", ncomp)
        W <- vector("list", ncomp)
        for (kk in 1:ncomp) {
            m[[kk]] <- (beta0 * m0 + Nk[kk] * xbar[[kk]]) / betas[kk]
            mlt <- (beta0 * Nk[kk]) / (beta0 + Nk[kk])
            W[[kk]] <- solve(W0inv + Nk[kk] * S[[kk]] + mlt * tcrossprod(xbar[[kk]] - m0))
        }
    }
    endTime <- proc.time()

    list(
        'alpha' = alpha,
        'beta' = betas,
        'mu' = m,
        'W' = W,  ## W is information matrix!
        'v' = v,
        'L' = L,
        'time' = endTime[3] - startTime[3]
    )
}
