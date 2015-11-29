#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "wig_funs.R",
             sep="/"))

## Parameters and Arguments ###################################################

spec <- matrix(c("input",       "i", 1, "character",
                 "output",      "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

## Do it ######################################################################

cover <- readBigWig(args$)

## Save as an RData ###########################################################

message("stroing output as RData")
save(cover, file=args[["output"]])

###############################################################################
