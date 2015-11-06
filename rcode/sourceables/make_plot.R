#!/usr/bin/Rscript

library(IRanges)

makePlot <- function(dyn, plot.start, plot.end)
{   # Given an already subsetted dynamic, and its start and end, plot it
    covs <- lapply(dyn,
                   function(x) coverage(x$originals))

    indel.covs <- lapply(dyn,
                         function(x) coverage(x$indels))

    initPlot(c(plot.start, plot.end),
             c(0, max(unlist(lapply(covs, as.vector)))))

    addCovs(covs[[1]], covs[[2]])
    addIndels(indel.covs[[1]], indel.covs[[2]])

    par.usr <- par("usr")
    par.ypc <- (par.usr[4] - par.usr[3]) / 100
    middle <- (par.usr[4] - par.usr[3]) * 0.5

    addLeft(dyn[[1]]$left.shifts,
            dyn[[2]]$left.shifts,
            middle,
            par.ypc)
    addRight(dyn[[1]]$right.shifts,
             dyn[[2]]$right.shifts,
             middle,
             par.ypc)
}

initPlot <- function(xlim, ylim)
    # Initiate an empty plot with the specified size
    plot(1,
         type='n',
         xlim=xlim,
         ylim=ylim,
         xlab="",
         ylab="",
         main="")

addCovs <- function(cov1, cov2)
{   # Plot coverages for the first and second experiment.
    # First experiment will be plotted as a grey area, and second experiment
    # wil be plotted as a black dotted line.
    lines(cov1, col="grey", type='h', lty=1, lwd=2)
    lines(cov2, col="black", lty=3, lwd=2)

    legend("topleft",
           c(paste("Ref1"), paste("Ref2")),
           col=c("grey", "black"),
           lty=c(1,3),
           lwd=c(3,3),
           bty="n",
           cex=0.8)
}

addIndels <- function(del.cov, ins.cov)
{   # Plot the coverages for insertions and delitons.
    # Insertions will be plotted as a green area and deletions as a red one.
    lines(del.cov, type='h', lwd=2, col="red")
    lines(ins.cov, type='h', lwd=2, col="green")

    legend("topleft",
           c("", "",
             "Inserted reads",
             "Removed reads"),
           col=c(NA, NA,
                 "green", "red"),
           lty=c(NA, NA, 1, 1),
           lwd=c(NA, NA, 3, 3),
           bty="n",
           cex=0.8)
}

dyadPos <- function(ran)
    # Find the center points of an IRanges
    round((start(ran) + end(ran)) / 2)


parseShifts <- function(x, y)
    # Helper function to compare pairs of shifts and put them in an
    # easily plotable format
    IRanges(start=dyadPos(x),
            end=dyadPos(x) + abs(dyadPos(y) - dyadPos(x)))

addLeft <- function(lefts1, lefts2, middle, par.ypc)
{   # Plot the upstream shifts as dark red arrows
    left.shifts <- parseShifts(lefts1, lefts2)
    db.left <- disjointBins(left.shifts)

    for (i in seq_along(left.shifts)) {
        x1 <- end(left.shifts)[i]
        y1 <- middle - par.ypc*db.left[i]
        x2 <- start(left.shifts)[i]
        y2 <- middle - par.ypc*db.left[i]

        arrows(x1, y1, x2, y2,
               length=0.05,
               col="darkred",
               angle=30,
               code=2,
               lwd=1)
    }

    legend("topleft",
           c("", "", "", "",
             "Upstream shift"),
           col=c(NA, NA, NA, NA,
                 "darkred"),
           lty=c(NA, NA, NA, NA, 1),
           lwd=c(NA, NA, NA, NA, 3),
           bty="n",
           cex=0.8)
}

addRight <- function(rights1, rights2, middle, par.ypc)
{   # Plot the downstream arrows as dark blue arrows
    right.shifts <- parseShifts(rights1, rights2)
    db.right <- disjointBins(right.shifts)

    for (i in seq_along(right.shifts)) {
        x1 <- start(right.shifts)[i]
        y1 <- middle + par.ypc*db.right[i]
        x2 <- end(right.shifts)[i]
        y2 <- middle + par.ypc*db.right[i]

        arrows(x1, y1, x2, y2,
               length=0.05,
               col="darkblue",
               angle=30,
               code=2,
               lwd=1)
    }

    legend("topleft",
           c("", "", "", "", "",
             "Downstream shift"),
           col=c(NA, NA, NA, NA, NA,
                 "darkblue"),
           lty=c(NA, NA, NA, NA, NA, 1),
           lwd=c(NA, NA, NA, NA, NA, 3),
           bty="n",
           cex=0.8)
}
