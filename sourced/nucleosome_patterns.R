SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourced"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

###############################################################################
# Top wrapper function ########################################################

nucleosomePatternsDF <- function(calls, cover=NULL, df, col.id="name",
                                 col.chrom="chrom", col.pos="pos",
                                 col.strand="strand", mc.cores=1, ...)
{
    n <- nrow(df)
    iterRows <- function(i, report.every=100) {
        if (i %% report.every == 0) {
            message(i, "/", n)
        }
        nucleosomePatterns(calls  = calls,
                           cover  = cover,
                           id     = df[i, col.id],
                           chrom  = df[i, col.chrom],
                           pos    = df[i, col.pos],
                           strand = df[i, col.strand])
    }
    do.call(rbind,
            xlapply(1:n,
                    iterRows,
                    report.every=10,
                    mc.cores=mc.cores))
}

# Actual main function ########################################################
nucleosomePatterns <- function(calls, cover=NULL, id, chrom, pos, strand="+",
                               window=300, p1.max.merge=3,
                               p1.max.downstream=20, open.thresh=215,
                               max.uncovered=150)
{   # This function is written in a monadic style. This allows to express in
    # a more compact way the fact that every computation might be the last one
    # to be perfomed is some condition is met.
    flow <- list(continue = TRUE,
                 state    = list(id     = id,
                                 chrom  = as.character(chrom),
                                 strand = as.character(strand),
                                 pos    = pos,
                                 p1.pos = NA,
                                 m1.pos = NA,
                                 dist   = NA,
                                 descr  = NA))

    fs <- rev(list(partial(lookCover, cover, window, max.uncovered),
                   partial(getNearby, calls, window),
                   partial(getP1, p1.max.downstream, p1.max.merge),
                   getM1,
                   partial(getDescr, open.thresh)))

    doIt <- do.call(compose,
                    lapply(fs,
                           procFun,
                           compose(bind,
                                   wrapFun)))

    makeDfRow(doIt(flow)$state)
}

###############################################################################
# Helper functions ############################################################

.mid <- function(x)
    floor((start(x)+end(x)) / 2)

.gimmeDist <- function(x, open.thresh)
    ifelse(is.na(x),        "-",
    ifelse(x > open.thresh, "open",
    ifelse(x < 120,         "overlap",
                            "close")))

makeDfRow <- function(ls)
{
    with(ls,
         data.frame(id     = id,
                    chrom  = chrom,
                    strand = strand,
                    pos    = pos,
                    p1.pos = p1.pos,
                    m1.pos = m1.pos,
                    dist   = dist,
                    descr  = descr))
}

detectNuc <- function (xs, pos, margin, strand, nucpos)
{
    a <- strand == "+" && nucpos == "p1"
    b <- strand == "-" && nucpos == "m1"
    c <- strand == "-" && nucpos == "p1"
    d <- strand == "+" && nucpos == "m1"

    after <- a || b
    before <- c || d

    if (after) {
        shift <- `-`
        comp <- `>`
        closest <- which.min
    } else if (before) {
        shift <- `+`
        comp <- `<`
        closest <- which.max
    }

    subxs <- xs[comp(.mid(xs), shift(pos, margin)), ]
    x <- subxs[closest(.mid(subxs)), ]
}

###############################################################################
# Functions to be used monadically

getNearby <- function (state, calls, window)
{
    mids <- .mid(calls)
    sel <- mids > (state$pos - window) &
           mids < (state$pos + window) &
           space(calls) == state$chrom
    if (any(sel)) {
        return(list(continue=TRUE,
                    vals=list(nearby=calls[sel, ])))
    } else {
        return(list(continue=TRUE,
                    vals=list()))
    }
}

getP1 <- function (state, p1.max.downstream, p1.max.merge)
{
    p1 <- detectNuc(state$nearby,
                    state$pos,
                    p1.max.downstream,
                    state$strand, "p1")

    if (!nrow(p1)) {
        return(list(continue=FALSE, vals=list(descr="+1_missing")))
    } else if (p1$nmerge > p1.max.merge) {
        return(list(continue=FALSE, vals=list(descr="+1_too_fuzzy")))
    } else {
        return(list(continue=TRUE,
                    vals=list(p1.pos   = .mid(p1),
                              p1.class = p1$class)))
    }
}

getM1 <- function (state)
{
    m1 <- detectNuc(state$nearby, state$p1.pos, 0, state$strand, "m1")

    if (!nrow(m1)) {
        return(list(continue=FALSE, vals=list(descr="-1_missing")))
    } else {
        return(list(continue=TRUE,
                    vals=list(m1.pos   = .mid(m1),
                              m1.class = m1$class)))
    }
}

getDescr <- function (res, open.thresh)
{
    dist <- abs(res$p1.pos - res$m1.pos)
    dist.class <- .gimmeDist(dist, open.thresh)
    descr <- paste(res$m1.class,
                   dist.class,
                   res$p1.class,
                   sep="-")
    return(list(continue=FALSE,
                vals=list(dist=dist,
                          descr=descr)))
}

lookCover <- function (state, cover, window, max.uncovered)
{
    outOfBounds <- function (pos, chrom, window, cover)
        (pos-window) < 0 | (pos+window) > length(cover[[chrom]])

    isUncovered <- function (pos, chrom, window, cover, max.uncovered)
        sum(cover[[chrom]][(pos-window):(pos+window)] == 0) > max.uncovered

    cover.problem <- `&&`(!is.null(cover),
                          `||`(outOfBounds(state$pos,
                                           state$chrom,
                                           window,
                                           cover),
                               isUncovered(state$pos,
                                           state$chrom,
                                           window,
                                           cover,
                                           max.uncovered)))
    if (cover.problem) {
        return(list(continue=FALSE, vals=list()))
    } else {
        return(list(continue=TRUE, vals=list()))
    }
}

###############################################################################
# Some declarations to make the monad work ####################################

bind <- function (f)
    # state -> flow -> flow -> flow
    function (flow)
        if (flow$continue) {
            return(f(flow$state))
        } else {
            return(flow)
        }

wrapFun <- function(f)
    # state -> vals -> state -> flow
    function(state)
        with(f(state),
             list(continue = continue,
                  state    = updateVals(state, vals)))

###############################################################################
