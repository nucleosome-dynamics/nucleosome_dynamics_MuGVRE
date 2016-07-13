#!/usr/bin/env Rscript

###############################################################################

library(IRanges)
library(plyr)

###############################################################################

isIn <- function (x, y)
    # Is x in y?
    start(x) <= end(y) & end(x) >= start(y)

myFilter <- function (x, f, ...)
    x[f(x, ...)]

selectReads <- function (reads, range.df)
    myFilter(reads,
             isIn,
             range(with(range.df,
                        IRanges(start,
                                end))))

###############################################################################

Gaussian <- function (k=1, m=0, sd=1)
{
    runGaussian <- function (x) {
        num <- (x - m) ** 2
        denum <- 2 * sd**2
        k * exp(-num/denum)
    }
    plotGaussian <- function (x, ...)
        plot(x, runGaussian(x), ...)
    vals <- list(k=k, m=m, sd=sd)
    toList <- function ()
        vals
    structure(c(vals,
                apply=runGaussian,
                plot=plotGaussian),
              class="Gaussian")
}

print.Gaussian <- function (x)
{
    ns <- c("k", "sd", "m")
    print(paste(mapply(paste,
                       ns,
                       lapply(lapply(ns,
                                     flip2args(`[[`),
                                     x),
                              round,
                              digits=3),
                       MoreArgs=list(sep="=")),
                collapse=" "))
}

as.list.Gaussian <- function (x)
    list(k=x$k,
         m=x$m,
         sd=x$sd)

as.data.frame.Gaussian <- function (x)
    as.data.frame(as.list(x))

###############################################################################

stripZeros <- function (x, a=0, tail=FALSE)
{   # remove consecutive zeroes from the beginning of a vector
    # (or from the end if tail=TRUE)
    if (tail) {
        x <- rev(x)
    }
    if (x[1] == a) {
        shift <- rle(x)$lengths[1]
    } else {
        shift <- 0
    }
    x <- x[-(1:shift)]
    if (tail) {
        x <- rev(x)
    }
    list(x=x, shift=shift)
}

fitIt <- function (f, tries=10, max.value=35000, max.sd=2000)
{
    isFitOk <- function (fit)
        fit$value <= max.value && abs(fit$par[2]) <= max.sd
    for (i in 1:tries) {
        par <- 10 * runif(3)
        fit <- optim(par,
                     f,
                     method="BFGS",
                     control=list(reltol=1e-9))
        ok <- isFitOk(fit)
        if (ok) {
            break
        }
    }
    if (ok) {
        return(fit$par)
    } else {
        return(NULL)
    }
}

#gaussFitCov <- function (xs)
#{
#    xs <- as.vector(xs)
#    a <- stripZeros(xs, 0)
#    xs <- a$x
#    shift <- a$shift
#    xs <- stripZeros(xs, 0, tail=TRUE)$x
#    f <- function (par) {
#        m <- par[1]
#        sd <- par[2]
#        k <- par[3]
#        rhat <- k * exp(-0.5 * ((seq_along(xs) - m)/sd)^2)
#        sum((xs - rhat)^2)
#    }
#    fit.par <- fitIt(f)
#    if (is.null(fit.par)) {
#        return(NULL)
#    } else {
#        return(do.call(Gaussian,
#                       list(m=fit.par[[1]] + shift,
#                            sd=abs(fit.par[[2]]),
#                            k=fit.par[[3]])))
#    }
#}

#gaussFitCov <- function (xs)
#{
#    xs <- as.vector(xs)
#    vals <- unlist(mapply(rep,
#                          seq_along(cov),
#                          xs,
#                          SIMPLIFY=FALSE))
#    Gaussian(m=mean(vals),
#             sd=sd(vals),
#             k=max(xs))
#}

gaussFitCov <- function (xs)
{
    counts <- as.numeric(xs)
    vals <- seq_along(counts)
    total <- sum(counts)

    m <- sum(counts * vals) / total
    sd <- sqrt(sum((m - vals)^2 * counts) / total)
    k <- max(counts)

    Gaussian(m=m, sd=sd, k=k)
}


###############################################################################

gaussFitChr <- function (calls, reads, ...)
    adply(calls,
          1,
          function (x)
              as.data.frame(gaussFitCov(coverage(selectReads(reads, x)))),
          ...)

doGaussFit <- function (calls, reads, ...)
{
    if (class(reads) == "RangedData") {
        rs <- ranges(reads)
        doChr <- function (chr.calls) {
            chr <- chr.calls[1, "seqname"]
            message(chr)
            chr.reads <- rs[[chr]]
            gaussFitChr(chr.calls, chr.reads, ...)
        }
        ddply(calls, "seqname", doChr)
    } else if (class(reads) == "IRanges") {
        gaussFitChr(calls, reads, ...)
    }
}

###############################################################################
