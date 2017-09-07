sampleFromMixMdl <- function (yvec, xmat, zmat, group, niter, nburn, nthin, id) {
    library(inline)
    library(Rcpp)
    library(rstan)

    gg <- ordered(as.character(group), levels = sort(unique(group)))
    group <- as.integer(gg)

    simList = list(
        nobs = length(yvec),
        nfixef = ncol(xmat),
        nranef = ncol(zmat),
        ngroup = length(unique(group)),
        xmat = xmat,
        zmat = zmat,
        group = group,
        yvec = yvec)

    seeds <- (1:5000) * as.numeric(gsub(":", "", substr(Sys.time(), 12, 19)))

    stanCode <- readChar("full_lme.stan", file.info("full_lme.stan")$size)
    startTime <- proc.time()
    mdl <- stan(model_code = stanCode, data = simList, iter = niter, warmup = nburn, chains = 1, thin = nthin,
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
