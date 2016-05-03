#!/usr/bin/env Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)
library(plyr)

where <- function () {
    spath <-parent.frame(2)$ofile

    if (is.null(spath)) {
        args <- commandArgs()
        filearg <- args[grep("^--file=", args)]
        fname <- strsplit(filearg, "=")[[1]][2]
    } else {
        fname <- spath
    }

    dirname(normalizePath(fname))
}

#SOURCE.DIR <- paste(where(), "../sourced", sep="/")
SOURCE.DIR <- "/home/rilla/nucleServ/sourced"
sourced <- c("gff_funs", "fp", "gauss_funs")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
     source(x)
}

## Command line arguments #####################################################

spec <- matrix(c("calls",  "a", 1, "character",
                 "reads",  "b", 1, "character",
                 "output", "c", 1, "character",
                 "start",  "d", 1, "integer",
                 "end",    "e", 1, "integer",
                 "chr",    "f", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)


args$calls="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
args$reads="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
args$start=1000
args$end=10000
args$chr="chrI"
#args$output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_gauss.gff"


## Subset #####################################################################

message("loading inputs")
calls <- readGff(args[["calls"]])
reads <- get(load(args[["reads"]]))

if (!is.null(args[["chr"]])) {
    message("subsetting data")

    calls <- subset(calls, seqname == args[["chr"]])
    reads <- ranges(reads)[[args[["chr"]]]]

    if (!is.null(args$start) && !is.null(args$end)) {
        calls <- subset(calls,
                        isIn(IRanges(start, end),
                             IRanges(args[["start"]],
                                     args[["end"]])))
        reads <- selectReads(reads, calls)
    }
}

message("perfroming the fittings")
gauss.df <- doGaussFit(calls, reads, .progress="text")
gauss.df$stiffess <- 1/gauss.df$sd

## Save output ################################################################

gauss.df$score_w <- NULL
gauss.df$score_h <- NULL
gauss.df$nmerge <- NULL

gff <- df2gff(gauss.df)
writeGff(gff, args[["output"]])

###############################################################################
