#!/usr/bin/Rscript

library(getopt)
library(NucDyn)

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

png(filename=args[["output"]],
    width=1000,
    height=500)
with(args,
     plotDynamics(get(load(input)),
                  chr=chr,
                  plot.range=c(start, end),
                  main=""))
dev.off()
