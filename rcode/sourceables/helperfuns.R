#!/usr/bin/Rscript

library(IRanges)
library(GenomicRanges)

source(paste(SOURCE.DIR,
             "fp.R",
             sep="/"))

# Miscelanious helper functions

vectorizedAll <- function(...)
    # Helper function that behaves as a vectorized version of the function
    # `all`
    Reduce(`&`, list(...))

getType <- function(score_w, score_h, thresh_w, thresh_h)
    # Given a vector of score w, a vector of score h, a threshold for w and a
    # threshold for h, return a vector that tells whether those scores are
    # classified as `well-positioned`, or `fuzzy`.
    ifelse(`&`(score_w > thresh_w,
               score_h > thresh_h),
           "W", "F")

checkInF <- function(f)
{   # Some checks on the input file specified
    if (is.null(f)) {
        stop("An input file must be specified")
    }
    if (!grepl("\\.bam$", f)) {
        stop("Input file must be in BAM format")
    }
    if (!file.exists(f)) {
        stop("Specified input file doesn't exist")
    }
}

checkType <- function(t)
{   # Some checks on the sequencing type specified
    if (is.null(t)) {
        stop("A type has to be specified")
    }
    if (!t %in% c("single", "paired")) {
        stop("Type must be either `single` or `paired`")
    }
}

checkOutF <- function(f)
{   # Some checks on the output file specified
    if (is.null(f)) {
        stop("An output file must be specified")
    }
}

makePlotable <- function(dyn)
{
    message("building structure to be saved for future plotting")

    useful.types <- c("originals", "right.shifts", "left.shifts", "indels")
    chrs <- unique(unlist(lapply(seqnames(set.a(dyn)), levels)))

    lapply(
        list(set.a(dyn), set.b(dyn)),
        function(set) {
            by.chrs <- lapply(
                chrs,
                function(chr)
                    as.list(ranges(set[seqnames(set) == chr, ])[useful.types])
            )
            names(by.chrs) <- chrs
            by.chrs
        }
    )
}

sortDfBy <- function(df, xs)
    do.call(compose,
            lapply(xs,
                   flip2args(partial),
                   orderBy))(df)

orderBy <- function(df, x)
    df[order(df[, x]), ]

getFirstTx <- function(x, df)
{
    entries <- subset(df, GENEID == x)
    f <- `[[`(list("+"=which.min,
                   "-"=which.max),
              unique(entries$TXSTRAND))
    entries[f(entries$TXSTART), ]
}

subMany <- function (patterns, replacements, x)
    do.call(compose,
            mapply(function(p, r) {force(p)
                                   force(r)
                                   function(x) sub(p, r, x)},
                   patterns,
                   replacements
    ))(x)

###############################################################################

.check.mc <- function (mc.cores) {
    lib <- "parallel"
    if (mc.cores > 1 && !lib %in% loadedNamespaces()) {
        warning("'",
                lib,
                "' library not available, switching to m.cores=1")
        return(1)
    } else {
        return(mc.cores)
    }
}

xlapply <- function(X, FUN, ..., mc.cores=1)
{   # Choose between lapply pr mclapply accordingly
    actual.cores <- .check.mc(mc.cores)
    if (actual.cores > 1) {
        mclapply(X=X,
                 FUN=FUN,
                 ...=...,
                 mc.cores=actual.cores)
    } else {
        lapply(X=X,
               FUN=FUN,
               ...=...)
    }
}

updateVals <- function (df, vals)
{
    for (i in names(vals)) {
        df[[i]] <- vals[[i]]
    }
    return(df)
}

dyadPos <- function (x)
    (start(x) + end(x))/2

iterDf <- function(df, fun, ...)
    lapply(1:nrow(df),
           function(i) do.call(fun,
                               c(unname(as.list(df[i, ])),
                                 list(...))))
