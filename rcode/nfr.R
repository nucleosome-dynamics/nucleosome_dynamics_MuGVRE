#!/usr/bin/env Rscript

# Given a nucleosome map, get the nucleosome-free regions, defined as linker
# ranges that fullfill some given conditions

# SOLIPSISTIC PUBLIC LICENSE
# Version 1, April 2013
# 
# Copyright (C) 2016 Ricard Illa Pujagut
# 
# Everyone is permitted to copy and distribute verbatim copies of
# this license document. Modified copies of this document are 
# permitted provided that they denounce BOTH the original AND their
# copy as mere sense data with no verifiable cause outside the mind.
# 
#                     SOLIPSISTIC PUBLIC LICENSE
#   TERMS AND CONDITIONS FOR COPYING, DISTRIBUTION AND MODIFICATION
# 
# 0. The term 'work' refers to the false sense-data distributed
#    with this license.
# 1. The term 'you' refers to the only being who verifiably exists.
# 2. The term 'author' refers to the set of delusions whereby you
#    incorrectly assign external agency to the work.
# 3. You may copy, modify and distribute the work without restrictions
#    provided that you do not believe the author exists, and provided
#    that you affirm publicly when referring to the work, or when
#    questioned or interrogated by beings who putatively exist, that
#    the work does not exist. 

## Imports ####################################################################

library(getopt)
library(IRanges)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
sourced <- c("helperfuns", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

# default parameters
defaults <- list(min.width = 400,
                 threshold = 110)

# parse arguments from the command line
spec <- matrix(c("input",     "i", 1, "character",
                 "output",    "o", 1, "character",
                 "minwidth",  "m", 1, "integer",
                 "threshold", "t", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("minwidth", "min.width", names(args))

# set all paremeters to be used
params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some declarations ##########################################################

# conditions that a linker has to satisfy for it to be considered a NFR.
# so far, its width has to be higher than the minimum width and lower than
# a threshold
conditions <- list(compose(partial(`>`, params[["threshold"]]), width),
                   compose(partial(`<`, params[["min.width"]]), width))

getInterRans <- function(r)
{   # Given an IRanges, return its negative ranges.
    # Use it to get the linker ranges of an IRanges representing a nucleosome
    # map
    redr <- reduce(r)
    IRanges(start=end(redr)[-length(redr)] + 1,
            end=start(redr)[-1] - 1)
}

## Read input #################################################################

message("-- reading input")
tab <- readGff(params[["input"]])  # read the nucleosome map as a gff3 input
names(tab) <- sub("seqname", "space", names(tab))
calls <- RangedData(tab)

## Do it ######################################################################

message("-- looking for Nucleosome Free regions")
nfr <- irLs2rd(lapply(
    ranges(calls),
    function(x) {
        linkers <- getInterRans(x)
        sel <- Reduce(`&`,
                      lapply(conditions,
                             application,
                             linkers))
        linkers[sel]
    }
))

## Store output ###############################################################

message("-- saving gff output")
# save the output as a gff3 too
writeGff(df2gff(rd2df(nfr),
                source="nucleR",
                feature="Nucleosme Free Region"),
         params[["output"]])

###############################################################################
