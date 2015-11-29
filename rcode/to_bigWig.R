#!/usr/bin/Rscriptq

## Imports ####################################################################

library(getopt)
library(IRanges)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "wig_funs.R",
             sep="/"))

## Parameters and Arguments ###################################################

spec <- matrix(c("input",  "i", 1, "character",
                 "output", "o", 1, "character",
                 "genome", "g", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

## Load RData #################################################################

cover <- get(load(args$input))

## Do it ######################################################################

writeBigWig(lapply(cover,
                   splitAtZeros),
            args$output,
            args$genome)

###############################################################################
