#!/usr/bin/Rscript

library(getopt)

spec <- matrix(c("input",  "i", 1, "character",
                 "output", "o", 1, "character",
                 "start",  "s", 1, "integer",
                 "end",    "e", 1, "integer",
                 "chr",    "c", 1, "character"),
               byrow=TRUE,
               ncol=4)

args <- getopt(spec)

if (any(!spec[, 1] %in% names(args))) {
    message("arguments:")
    for (x in spec[, 1]) {
        message("\t--", x)
    }
    q("no")
}

SOURCE.DIR <- "/home/rilla/nucleServ/sourced"
source(paste(SOURCE.DIR,
             "plot_subset.R",
             sep="/"))
source(paste(SOURCE.DIR,
             "make_plot.R",
             sep="/"))

library(IRanges)

dyn <- get(load(args$input))
subdyn <- subsetDyn(dyn, args$chr, args$start, args$end)

png(filename=args[["output"]], width=1000, height=500)
makePlot(subdyn, args$start, args$end)
dev.off()
