#!/usr/bin/env Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
sourced <- c("helperfuns", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

defaults <- list(min.width = 400,
                 threshold = 110)

spec <- matrix(c("input",     "i", 1, "character",
                 "output",    "o", 1, "character",
                 "minwidth",  "m", 1, "integer",
                 "threshold", "t", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("minwidth", "min.width", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some declarations ##########################################################

conditions <- list(compose(partial(`>`, params[["threshold"]]), width),
                   compose(partial(`<`, params[["min.width"]]), width))

getInterRans <- function(r)
{
    redr <- reduce(r)
    IRanges(start=end(redr)[-length(redr)] + 1,
            end=start(redr)[-1] - 1)
}

## Read input #################################################################

message("-- reading input")
tab <- readGff(params[["input"]])
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
writeGff(df2gff(rd2df(nfr),
                source="nucleR",
                feature="Nucleosme Free Region"),
         params[["output"]])

###############################################################################
