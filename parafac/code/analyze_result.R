rm(list = ls())
setwd("/Shared/ssrivastva/wasp/parafac/result/")

library(KernSmooth)
library(matrixStats)
library(RColorBrewer)
colors <- brewer.pal(6, "Set1")

full <- list()
for (cc in 1:10) {
    full[[cc]] <- vector("list", 20)
    names(full[[cc]]) <- paste0("dim", 1:20)
    for (dd in 1:20) {
        dat <- read.table(paste0("full/res_cv_", cc, "_dim_", dd, ".csv"), sep = ",", header = FALSE)
        colnames(dat) <- c("p1", "p2")
        full[[cc]][[dd]] <- dat
    }
}

sub10 <- list()
for (cc in 1:10) {
    sub10[[cc]] <- vector("list", 20)
    names(sub10[[cc]]) <- paste0("dim", 1:20)
    for (dd in 1:20) {
        sub10[[cc]][[dd]] <- vector("list", 10)
        names(sub10[[cc]][[dd]]) <- paste0("sub", 1:10)
        for (kk in 1:10) {
            dat <- read.table(paste0("sub10/samp/csv/res_cv_", cc, "_sub_", kk, "_dim_", dd, "_k10.csv"), sep = ",", header = FALSE)
            colnames(dat) <- c("p1", "p2")
            sub10[[cc]][[dd]][[kk]] <- dat
        }
    }
}

sub5 <- list()
for (cc in 1:10) {
    sub5[[cc]] <- vector("list", 20)
    names(sub5[[cc]]) <- paste0("dim", 1:20)
    for (dd in 1:20) {
        sub5[[cc]][[dd]] <- vector("list", 5)
        names(sub5[[cc]][[dd]]) <- paste0("sub", 1:5)
        for (kk in 1:5) {
            dat <- read.table(paste0("sub5/samp/csv/res_cv_", cc, "_sub_", kk, "_dim_", dd, "_k5.csv"), sep = ",", header = FALSE)
            colnames(dat) <- c("p1", "p2")
            sub5[[cc]][[dd]][[kk]] <- dat
        }
    }
}

wasp10 <- list()
for (cc in 1:10) {
    wasp10[[cc]] <- list()
    for (dd in 1:20) {
        vec1 <- rowMeans(do.call(cbind, lapply(sub10[[cc]][[dd]], function(x) quantile(x[ , "p1"], probs = seq(0, 1, length = 10000)))))
        vec2 <- rowMeans(do.call(cbind, lapply(sub10[[cc]][[dd]], function(x) quantile(x[ , "p2"], probs = seq(0, 1, length = 10000)))))
        wasp10[[cc]][[dd]] <- cbind("p1" = vec1, "p2" = vec2)
    }
}

wasp5 <- list()
for (cc in 1:10) {
    wasp5[[cc]] <- list()
    for (dd in 1:20) {
        vec1 <- rowMeans(do.call(cbind, lapply(sub5[[cc]][[dd]], function(x) quantile(x[ , "p1"], probs = seq(0, 1, length = 10000)))))
        vec2 <- rowMeans(do.call(cbind, lapply(sub5[[cc]][[dd]], function(x) quantile(x[ , "p2"], probs = seq(0, 1, length = 10000)))))
        wasp5[[cc]][[dd]] <- cbind("p1" = vec1, "p2" = vec2)
    }
}

xing5 <- list()
xing10 <- list()
cons5 <- list()
cons10 <- list()
for (cc in 1:10) {
    xdat <- readRDS(paste0("comp/sub5/marg/xing_cv_", cc, "_k5.rds"))
    cdat <- readRDS(paste0("comp/sub5/marg/cons_cv_", cc, "_k5.rds"))
    xing5[[cc]] <- xdat$marg
    cons5[[cc]] <- cdat$marg
    xdat <- readRDS(paste0("comp/sub10/marg/xing_cv_", cc, "_k10.rds"))
    cdat <- readRDS(paste0("comp/sub10/marg/cons_cv_", cc, "_k10.rds"))
    xing10[[cc]] <- xdat$marg
    cons10[[cc]] <- cdat$marg
}


library(KernSmooth)
library(matrixStats)

res <- vector("list", 7)
names(res) <- c("full", "wasp5", "wasp10", "xing5", "xing10", "cons5", "cons10")

for (cc in 1:10) {
    res[["full"]][[cc]] <- list()
    res[["wasp5"]][[cc]] <- list()
    res[["wasp10"]][[cc]] <- list()
    res[["xing5"]][[cc]] <- list()
    res[["xing10"]][[cc]] <- list()
    res[["cons5"]][[cc]] <- list()
    res[["cons10"]][[cc]] <- list()
    for (dd in 1:20) {
        res[["full"]][[cc]][[dd]] <- list()
        res[["wasp5"]][[cc]][[dd]] <- list()
        res[["wasp10"]][[cc]][[dd]] <- list()
        res[["xing5"]][[cc]][[dd]] <- list()
        res[["xing10"]][[cc]][[dd]] <- list()
        res[["cons5"]][[cc]][[dd]] <- list()
        res[["cons10"]][[cc]][[dd]] <- list()
        for (pp in 1:2) {
            if (pp == 1) {
                rr <- range(c(full[[cc]][[dd]][ , pp], wasp5[[cc]][[dd]][ , pp], wasp10[[cc]][[dd]][ , pp],
                              xing5[[cc]][[dd]], xing10[[cc]][[dd]], cons5[[cc]][[dd]], cons10[[cc]][[dd]]))
            } else {
                rr <- range(c(full[[cc]][[dd]][ , pp], wasp5[[cc]][[dd]][ , pp], wasp10[[cc]][[dd]][ , pp],
                              1 - xing5[[cc]][[dd]], 1 - xing10[[cc]][[dd]],
                              1 - cons5[[cc]][[dd]], 1 - cons10[[cc]][[dd]]))
            }
            fbw <- dpik(full[[cc]][[dd]][ , pp], range.x = rr)
            wbw1 <- dpik(wasp5[[cc]][[dd]][ , pp], range.x = rr)
            wbw2 <- dpik(wasp10[[cc]][[dd]][ , pp], range.x = rr)
            dens1 <- bkde(full[[cc]][[dd]][ , pp], bandwidth = fbw, range.x = rr)
            dens2 <- bkde(wasp5[[cc]][[dd]][ , pp], bandwidth = wbw1, range.x = rr)
            dens3 <- bkde(wasp10[[cc]][[dd]][ , pp], bandwidth = wbw2, range.x = rr)
            if (pp == 1) {
                xbw1 <- dpik(xing5[[cc]][[dd]], range.x = rr)
                xbw2 <- dpik(xing10[[cc]][[dd]], range.x = rr)
                cbw1 <- dpik(cons5[[cc]][[dd]], range.x = rr)
                cbw2 <- dpik(cons10[[cc]][[dd]], range.x = rr)
                dens4 <- bkde(xing5[[cc]][[dd]], bandwidth = xbw1, range.x = rr)
                dens5 <- bkde(xing10[[cc]][[dd]], bandwidth = xbw2, range.x = rr)
                dens6 <- bkde(cons5[[cc]][[dd]], bandwidth = cbw1, range.x = rr)
                dens7 <- bkde(cons10[[cc]][[dd]], bandwidth = cbw2, range.x = rr)
            } else {
                xbw1 <- dpik(1 - xing5[[cc]][[dd]], range.x = rr)
                xbw2 <- dpik(1 - xing10[[cc]][[dd]], range.x = rr)
                cbw1 <- dpik(1 - cons5[[cc]][[dd]], range.x = rr)
                cbw2 <- dpik(1 - cons10[[cc]][[dd]], range.x = rr)
                dens4 <- bkde(xing5[[cc]][[dd]], bandwidth = xbw1, range.x = rr)
                dens5 <- bkde(xing10[[cc]][[dd]], bandwidth = xbw2, range.x = rr)
                dens6 <- bkde(cons5[[cc]][[dd]], bandwidth = cbw1, range.x = rr)
                dens7 <- bkde(cons10[[cc]][[dd]], bandwidth = cbw2, range.x = rr)
            }
            res[["full"]][[cc]][[dd]][[pp]] <- dens1
            res[["wasp5"]][[cc]][[dd]][[pp]] <- dens2
            res[["wasp10"]][[cc]][[dd]][[pp]] <- dens3
            res[["xing5"]][[cc]][[dd]][[pp]] <- dens4
            res[["xing10"]][[cc]][[dd]][[pp]] <- dens5
            res[["cons5"]][[cc]][[dd]][[pp]] <- dens6
            res[["cons10"]][[cc]][[dd]][[pp]] <- dens7
        }
    }
}

acc <- array(NA,
             dim = c(10, 20, 6),
             dimnames = list(paste0("cv", 1:10),
                             paste0("dim", 1:20),
                             c("wasp5", "wasp10","xing5", "xing10", "cons5", "cons10")
                             ))

for (cc in 1:10) {
    for (dd in 1:20) {
        acc[cc, dd, 1] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["wasp5"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
        acc[cc, dd, 2] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["wasp10"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
        acc[cc, dd, 3] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["xing5"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
        acc[cc, dd, 4] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["xing10"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
        acc[cc, dd, 5] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["cons5"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
        acc[cc, dd, 6] <- (1 - sum(abs(res[["full"]][[cc]][[dd]][[pp]]$y  - res[["cons10"]][[cc]][[dd]][[pp]]$y) * diff(res[["full"]][[cc]][[dd]][[pp]]$x)[1]) / 2)
    }
}

tblMn <- apply(acc[, c(2, 4, 12, 14), ], 2:3, mean)
tblSd <- apply(acc[, c(2, 4, 12, 14), ], 2:3, sd)

accTbl <- matrix(paste0(format(round(tblMn, 2), nsmall = 2), " (", format(round(tblSd, 2), nsmall = 2), ")"), nrow = 4)
rownames(accTbl) <- rownames(tblMn)
colnames(accTbl) <- colnames(tblMn)

xtable::xtable(accTbl[ , c(rev(c(1, 3, 5)), rev(c(2, 4, 6)))])

### plots

#### density

rnames <- c(expression("pr("~x[2] == 1~")"), expression("pr("~x[4] == 1~")"), expression("pr("~x[12] == 1~")"), expression("pr("~x[14] == 1~")"))

pdf("~/wasp/parafac/result/img/marg_para.pdf", 40, 12)
par(mfrow = c(2, 4))
par(cex = 1)
par(mar = c(2, 4, 0, 2), oma = c(1, 1, 0.4, 0.4))
par(tcl = -0.02)
par(mgp = c(2, 0.6, 0))
cc <- 10; pp <- 1;
cnt <- 0
for (dd in c(2, 4, 12, 14)) {
    cnt <- cnt + 1
    rr <- range(c(res[["full"]][[cc]][[dd]][[pp]]$y, res[["wasp5"]][[cc]][[dd]][[pp]]$y,  res[["wasp10"]][[cc]][[dd]][[pp]]$y, res[["xing5"]][[cc]][[dd]][[pp]]$y, res[["xing10"]][[cc]][[dd]][[pp]]$y, res[["cons5"]][[cc]][[dd]][[pp]]$y, res[["cons10"]][[cc]][[dd]][[pp]]$y))
    ## k = 5
    plot(res[["full"]][[cc]][[dd]][[pp]]$x, res[["full"]][[cc]][[dd]][[pp]]$y, lty = "solid", ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr, lwd = 8, col = colors[1])
    xxlab <- axTicks(1)
    yylab <- axTicks(2)
    axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 2.5, las = 2)
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(format(round(xxlab, 4)), at = round(xxlab, 4), side = 1, line = 1, cex = 2.5)
    grid(lwd = 5)
    box(col = "grey40", lwd = 5)
    lines(res[["cons5"]][[cc]][[dd]][[pp]]$x, res[["cons5"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[2])
    lines(res[["xing5"]][[cc]][[dd]][[pp]]$x, res[["xing5"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[3])
    lines(res[["wasp5"]][[cc]][[dd]][[pp]]$x, res[["wasp5"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[5])
    mtext(rnames[cnt], side = 3, line = -4.5, adj = 0.1, cex = 4)
    mtext("k = 5", side = 3, line = -7.5, adj = 0.1, cex = 4)
    if (dd == 2) {
        legend("topright", c("MCMC", "CMC", "SDP", "WASP"), lty = c(rep("solid", 4)), col = c(colors[1:3], colors[5]),
               lwd = 8, bty = "n", cex = 2.5)
    }
}
cnt <- 0
for (dd in c(2, 4, 12, 14)) {
    cnt <- cnt + 1
    rr <- range(c(res[["full"]][[cc]][[dd]][[pp]]$y, res[["wasp5"]][[cc]][[dd]][[pp]]$y,  res[["wasp10"]][[cc]][[dd]][[pp]]$y, res[["xing5"]][[cc]][[dd]][[pp]]$y, res[["xing10"]][[cc]][[dd]][[pp]]$y, res[["cons5"]][[cc]][[dd]][[pp]]$y, res[["cons10"]][[cc]][[dd]][[pp]]$y))
    ## k = 10
    plot(res[["full"]][[cc]][[dd]][[pp]]$x, res[["full"]][[cc]][[dd]][[pp]]$y, lty = "solid", ylab = NA, xlab = NA, axes = FALSE, type = "l", ylim = rr, lwd = 8, col = colors[1])
    xxlab <- axTicks(1)
    yylab <- axTicks(2)
    axis(side = 2, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(format(yylab), at = yylab, side = 2, line = 0.5, cex = 2.5, las = 2)
    axis(side = 1, tck = -0.01, lwd = 3, cex = 2, labels = NA)
    mtext(format(round(xxlab, 4)), at = round(xxlab, 4), side = 1, line = 1, cex = 2.5)
    grid(lwd = 5)
    box(col = "grey40", lwd = 5)
    lines(res[["cons10"]][[cc]][[dd]][[pp]]$x, res[["cons10"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[2])
    lines(res[["xing10"]][[cc]][[dd]][[pp]]$x, res[["xing10"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[3])
    lines(res[["wasp10"]][[cc]][[dd]][[pp]]$x, res[["wasp10"]][[cc]][[dd]][[pp]]$y, lty = "solid", lwd = 8, col = colors[5])
    mtext(rnames[cnt], side = 3, line = -4.5, adj = 0.1, cex = 4)
    mtext("k = 10", side = 3, line = -7.5, adj = 0.1, cex = 4)
    if (dd == 2) {
        legend("topright", c("MCMC", "CMC", "SDP", "WASP"), lty = c(rep("solid", 4)), col = c(colors[1:3], colors[5]),
               lwd = 8, bty = "n", cex = 2.5)
    }
}
dev.off()

#### time

fullTime <- list()
for (cc in 1:10) {
    fullTime[[cc]] <- read.table(paste0("full/time_", cc, ".csv"), sep = ",", header = FALSE)
}

subTime10 <- list()
for (cc in 1:10) {
    subTime10[[cc]] <- vector("list", 10)
    for (kk in 1:10) {
        subTime10[[cc]][[kk]] <- read.table(paste0("sub10/samp/time_cv_", cc, "_sub_", kk, "_k10.csv"), sep = ",", header = FALSE)
    }
}

subTime5 <- list()
for (cc in 1:10) {
    subTime5[[cc]] <- vector("list", 5)
    for (kk in 1:5) {
        subTime5[[cc]][[kk]] <- read.table(paste0("sub5/samp/time_cv_", cc, "_sub_", kk, "_k5.csv"), sep = ",", header = FALSE)
    }
}

compTime10 <- list()
for (cc in 1:10) {
    compTime10[[cc]] <- vector("list", 10)
    for (kk in 1:10) {
        compTime10[[cc]][[kk]] <- read.table(paste0("../result/comp/sub10/samp/time_cv_", cc, "_sub_", kk, "_k10.csv"), sep = ",", header = FALSE)
    }
}

compTime5 <- list()
for (cc in 1:10) {
    compTime5[[cc]] <- vector("list", 5)
    for (kk in 1:5) {
        compTime5[[cc]][[kk]] <- read.table(paste0("../result/comp/sub5/samp/time_cv_", cc, "_sub_", kk, "_k5.csv"), sep = ",", header = FALSE)
    }
}

xingTime5 <- list()
consTime5 <- list()
for (cc in 1:10) {
    xingTime5[[cc]] <- readRDS(paste0("comp/sub5/marg/xing_cv_", cc, "_k5.rds"))$time
    consTime5[[cc]] <- readRDS(paste0("comp/sub5/marg/cons_cv_", cc, "_k5.rds"))$time
}

xingTime10 <- list()
consTime10 <- list()
for (cc in 1:10) {
    xingTime10[[cc]] <- readRDS(paste0("comp/sub10/marg/xing_cv_", cc, "_k10.rds"))$time
    consTime10[[cc]] <- readRDS(paste0("comp/sub10/marg/cons_cv_", cc, "_k10.rds"))$time
}

rtime <- list(log10(unlist(fullTime)),
              log10(sapply(compTime5, function(x) mean(unlist(x))) + sapply(consTime5, mean)),
              log10(sapply(compTime10, function(x) mean(unlist(x))) + sapply(consTime10, mean)),
              log10(sapply(compTime5, function(x) mean(unlist(x))) + sapply(xingTime5, mean)),
              log10(sapply(compTime10, function(x) mean(unlist(x))) + sapply(xingTime10, mean)),
              log10(unlist(lapply(subTime5, function(x) mean(unlist(x))))),
              log10(unlist(lapply(subTime10, function(x) mean(unlist(x)))))
              )

pdf("~/wasp/parafac/result/img/para_time.pdf", 6, 8)
par(cex = 1)
par(mar = c(6.7, 0, 0, 0), oma = c(3, 5, 0.2, 0.2))
par(tcl = -0.1)
par(mgp = c(2, 0.6, 0))

boxplot(rtime, ylab = NA, axes = FALSE, lwd = 2, ylim = c(4, 6),
        boxlwd = 4, boxwex = 0.4, whisklty = 1, whisklwd = 3, staplelty = 1,
        staplelwd = 3, medlwd = 4, outcex = 1.5)
grid(lwd=3)
box(col = "grey40", lwd = 4)
axis(side = 1, tck = -.01, labels = NA)
mtext(c("MCMC", "CMC (k=5)", "CMC (k=10)", "SDP (k=5)",
        "SDP (k=10)", "WASP (k=5)", "WASP (k=10)"),
      at = c(1:7), side = 1, line = 0, cex = 1.9, las = 2)
axis(side = 2, tck = -0.01, labels = NA, lwd = 3)
mtext(format(seq(4, 6, by = 0.5)), at = seq(4, 6, by = 0.5), side = 2, line = 0.3, las = 1, cex = 2)
mtext(expression(log[10] * " Seconds"), side = 2, outer = TRUE, cex = 2.5, line = 2.5)
dev.off()
