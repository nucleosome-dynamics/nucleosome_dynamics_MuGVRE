#!/usr/bin/Rscript

library(ggplot2)
library(IRanges)
library(plotly)
library(plyr)
library(reshape2)

dyadPos <- function(ran)
    # Find the center points of an IRanges
    round((start(ran) + end(ran)) / 2)

makeEqual <- function (xs) {
    maxlen <- max(sapply(xs, length))
    lapply(
        xs,
        function (x)
            c(x, rep(0, maxlen - length(x)))
    )
}

buildArrowDf <- function (xs, ys, middle, par.ypc, type="left") {
    if (length(xs) == 0 | length(ys) == 0) {
        data.frame(x0=vector(),
                   x1=vector(),
                   y0=vector(),
                   y1=vector(),
                   variable=factor(c(), levels=type))
    } else {
        start <- dyadPos(xs)
        end <- dyadPos(ys)
        db <- disjointBins(IRanges(start = mapply(min, start, end),
                                   end   = mapply(max, start, end)))
        if (type == "left") {
            vpos <- middle - par.ypc*db
        } else if (type == "right") {
            vpos <- middle + par.ypc*db
        }
        data.frame(x0=start,
                   x1=end,
                   y0=vpos,
                   y1=vpos,
                   variable=type)
    }
}

arrowHead <- function (x0, y0, l, par.ypc, forward=FALSE) {
    a <- sin(pi/4) * l
    if (forward) {
        list(list(x0=x0, y0=y0, x1=x0 - 10, y1=y0 + a*par.ypc),
             list(x0=x0, y0=y0, x1=x0 - 10, y1=y0 - a*par.ypc))
    } else {
        list(list(x0=x0, y0=y0, x1=x0 + 10, y1=y0 + a*par.ypc),
             list(x0=x0, y0=y0, x1=x0 + 10, y1=y0 - a*par.ypc))
    }
}

addArrows <- function (sdf, par.ypc, forward=TRUE) {
    if (nrow(sdf) > 0) {
        adply(
            sdf,
            1,
            function (row) {
                head <- arrowHead(row$x1,
                                  row$y1,
                                  0.5,
                                  par.ypc,
                                  forward=forward)
                var <- as.vector(row$variable)
                head.rows <- lapply(lapply(head,
                                           c,
                                           variable=var),
                                    as.data.frame)
                do.call(rbind, c(list(row), head.rows))
            }
        )
    } else {
        data.frame(x0=0,
                   y0=0,
                   x1=0,
                   y1=0,
                   variable=factor(levels(lsdf$variable)))
    }
}

buildGgplot <- function (mdf, shdf, plot.start, plot.end) {
    ggplot(mdf, aes(x=x, y=value)) +
        xlim(plot.start, plot.end) +
        xlab("") +
        ylab("") +
        geom_area(aes(fill=variable,
                      color=variable,
                      linetype=variable),
                  position="identity") +
        geom_segment(data=shdf,
                     mapping=aes(x=x0, y=y0, xend=x1, yend=y1,
                                 color=variable)) +
        scale_color_manual(name="",
                           values=c("setA"="gray",
                                    "setB"="black",
                                    "ins"="#00FF0090",
                                    "dels"="#FF000090",
                                    "left"="darkred",
                                    "right"="darkblue")) +
        scale_fill_manual(name="",
                          values=c("setA"="#C0C0C0",
                                   "setB"="#00000000",
                                   "ins"="#00FF0090",
                                   "dels"="#FF000095",
                                   "left"="darkred",
                                   "right"="darkblue")) +
        scale_linetype_manual(name="",
                              values=c("setA"="solid",
                                       "setB"="dashed",
                                       "ins"="solid",
                                       "dels"="solid",
                                       "left"="solid",
                                       "right"="solid")) +
        theme_bw()
}

ggplot2widget <- function (p) {
    pl <- plotly_build(p)
    file.remove("Rplots.pdf")
    pl[["data"]][[1]][["name"]] <- "Coverage 1"
    pl[["data"]][[2]][["name"]] <- "Coverage 2"
    pl[["data"]][[3]][["name"]] <- "Deletions"
    pl[["data"]][[4]][["name"]] <- "Insertions"
    pl[["data"]][[5]][["name"]] <- "Upstream shifts"
    pl[["data"]][[6]][["name"]] <- "Downstream shifts"
    as.widget(pl)
}

makeShifts <- function (dyn, middle, par.ypc) {
    lsdf <- buildArrowDf(dyn[[1]]$left.shifts,
                         dyn[[2]]$left.shifts,
                         middle,
                         par.ypc,
                         "left")
    rsdf <- buildArrowDf(dyn[[1]]$right.shifts,
                         dyn[[2]]$right.shifts,
                         middle,
                         par.ypc,
                         "right")
    lsdf <- addArrows(lsdf, par.ypc, forward=FALSE)
    rsdf <- addArrows(rsdf, par.ypc, forward=TRUE)
    rbind(lsdf, rsdf)
}

makeCovs <- function (dyn) {
    cols <- list(setA=dyn[[1]]$originals,
                 setB=dyn[[2]]$originals,
                 dels=dyn[[1]]$indels,
                 ins=dyn[[2]]$indels)
    covs <- lapply(cols, function(x) as.vector(coverage(x)))
    valsdf <- as.data.frame(makeEqual(covs))
    valsdf$x <- 1:nrow(valsdf)
    covs <- lapply(cols, function(x) as.vector(coverage(x)))
    valsdf <- as.data.frame(makeEqual(covs))
    valsdf$x <- 1:nrow(valsdf)
    melt(valsdf, id.vars="x")
}

buildPlot <- function (dyn, start, end) {
    mdf <- makeCovs(dyn)
    maxy <- max(mdf$value)
    middle <- maxy / 2
    par.ypc <- maxy / 100
    shdf <- makeShifts(subdyn, middle, par.ypc)
    gg <- buildGgplot(mdf, shdf, start, end)
    ggplot2widget(gg)
}
