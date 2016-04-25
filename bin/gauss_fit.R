#!/usr/bin/env Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)
library(plyr)

SOURCE.DIR <- "/home/rilla/nucleServ/sourced"
sourced <- c("gff_funs", "fp", "gauss_funs")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
     source(x)
}

## Command line arguments #####################################################

spec <- matrix(c("calls", "a", 1, "character",
                 "reads", "b", 1, "character",
                 "start", "c", 1, "integer",
                 "end",   "d", 1, "integer",
                 "chr",   "e", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

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

## Save output ################################################################

gauss.df$score_w <- NULL
gauss.df$score_h <- NULL
gauss.df$nmerge <- NULL

gff <- df2gff(gauss.df)
writeGff(gff, args[["output"]])

###############################################################################
