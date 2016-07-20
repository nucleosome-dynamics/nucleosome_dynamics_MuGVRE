#!/usr/bin/Rscript

###############################################################################

library(IRanges)

source(paste(SOURCE.DIR,
             "fp.R",
             sep="/"))

###############################################################################

simpleMerger <- function (x) {
    readsInvolved <- sum(x$readsInvolved)
    nreads <- sum(x$nreads)
    totalReads <- x[1, "totalReads"]
    data.frame(coord = round(mean(x$coord)),
               start = min(x$start),
               end   = max(x$end),
               type = x[1, "type"],
               nuc  = x[1, "nuc"],
               chr  = x[1, "chr"],
               nreads        = nreads,
               totalReads    = totalReads,
               readsInvolved = readsInvolved,
               freads        = nreads / totalReads,
               hreads        = nreads / readsInvolved)
}

mergeKind <- function (x, nuc.width) {
    ovlp <- findOverlaps(IRanges(start=x$coord - round(nuc.width/2),
                                 width=nuc.width))
    i <- queryHits(ovlp)[queryHits(ovlp) != subjectHits(ovlp)]
    if (length(i) > 1) {
        y <- x[i, ]
        red <- reduce(IRanges(start=y$coord - round(nuc.width/2),
                              width=nuc.width))
        merged <- do.call(rbind,
                          mapply(function (s, e)
                                     ddply(y[y$start >= s & y$end <= e, ],
                                           "nuc",
                                           simpleMerger),
                                 start(red),
                                 end(red),
                                 SIMPLIFY=FALSE))
        rbind(x[-i, ], merged)
    } else {
        x
    }
}

comb2shifts <- function (r, l) {
    coord <- mean(r$coord, l$coord)
    nreads <- r$nreads + l$nreads
    start <- min(r$start, l$start)
    end <- max(r$end, l$end)

    if (r$coord < l$coord) {
        type <- "decrease in fuzziness"
    } else if (r$coord > l$coord) {
        type <- "increase in fuzziness"
    } else {
        type <- NA
    }

    nuc <- r$nuc
    totalReads <- r$totalReads
    readsInvolved <- r$readsInvolved + l$readsInvolved
    totalReads <- r$totalReads
    chr <- r$chr
    freads <- nreads / totalReads
    hreads <- nreads / readsInvolved

    data.frame(coord         = coord,
               nreads        = nreads,
               start         = start,
               end           = end,
               type          = type,
               nuc           = nuc,
               totalReads    = totalReads,
               freads        = freads,
               readsInvolved = readsInvolved,
               hreads        = hreads,
               chr           = chr)
}

shiftCombiner <- function (i, j, rsh, lsh, same.magnitude) {
    r <- rsh[i, ]
    l <- lsh[j, ]
    nrs <- c(r$nreads, l$nreads)
    magnitude <- max(nrs) / min(nrs) <= params$same.magnitude
    magnitude <- ifelse(is.na(magnitude), FALSE, magnitude)
    nuc <- r$nuc == l$nuc && r$nuc != 0
    if (magnitude && nuc) {  # we'll combine them
        comb2shifts(r, l)
    } else {
        rbind(r, l)
    }
}

shiftChrMerger <- function (df, nuc.width, same.magnitude) {
    ii <- df$type == "SHIFT +"
    jj <- df$type == "SHIFT -"

    rsh <- df[ii, ]
    lsh <- df[jj, ]

    rran <- IRanges(start=rsh$start - round(nuc.width/2), width=nuc.width)
    lran <- IRanges(start=lsh$start - round(nuc.width/2), width=nuc.width)
    ovlp <- findOverlaps(rran, lran)

    merged <- do.call(rbind,
                      mapply(shiftCombiner,
                             queryHits(ovlp),
                             subjectHits(ovlp),
                             MoreArgs=list(rsh, lsh, same.magnitude),
                             SIMPLIFY=FALSE))

    rbind(df[!(ii | jj), ], merged)
}

hsSorter <- function (df) {
    df <- df[order(df$coord), ]
    df <- df[order(df$chr), ]
    unrowname(df)
}

###############################################################################

newCombiner <- function (hs, nuc.width, same.magnitude) {
    f <- compose(hsSorter,
                 partial(ddply,
                         .(chr),
                         shiftChrMerger,
                         nuc.width,
                         same.magnitude,
                         .progress="text"),
                 partial(ddply,
                         .(type, chr),
                         mergeKind,
                         nuc.width,
                         .progress="text"),
                 function (x) x[-grep("^CONTAINED ", x$type), ])
    f(hs)
}

###############################################################################
