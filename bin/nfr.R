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

SOURCE.DIR <- paste(where(), "../sourced", sep="/")
sourced <- c("helperfuns", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

# default parameters
defaults <- list(min.width = 400,
                 threshold = 110)

# parse arguments from the command line
spec <- matrix(c("input",      "i", 1, "character",
                 "output",     "o", 1, "character",
                 "minwidth",   "m", 1, "integer",
                 "threshold",  "t", 1, "integer",
                 "readsInput", "r", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("minwidth", "min.width", names(args))

# set all paremeters to be used
params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

params$reads <- "/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep2_00m_G1.bam.RData"
params$input <- "/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/run041/NR_rep2_00m_G1.gff"

## Some declarations ##########################################################

# conditions that a linker has to satisfy for it to be considered a NFR.
# so far, its width has to be higher than the minimum width and lower than
# a threshold
conditions <- list(
    # compose(partial(`>`, params[["threshold"]]), width),
    # compose(partial(`<`, params[["min.width"]]), width)
    function (x) width(x) > params[["min.width"]]
)

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
nfr <- rd2df(irLs2rd(lapply(
    ranges(calls),
    function(x) {
        linkers <- getInterRans(x)
        sel <- Reduce(`&`,
                      lapply(conditions,
                             application,
                             linkers))
        linkers[sel]
    }
)))

rmUnderCovered <- function (df, reads, thresh) {

    getCovSum <- function (from, to, cov) {
        sub.cov <- cov[from:to]
        sum(sub.cov == 0) / length(sub.cov)
    }

    doByChr <- function (df, covs, thresh) {
        chr <- as.vector(df[1, "seqname"])
        chr.cov <- covs[[chr]]
        covs.ratio <- mapply(getCovSum,
                             df$start,
                             df$end,
                             MoreArgs=list(chr.cov))
        df[covs.ratio < thresh, ]
    }

    covs <- lapply(coverage(reads), as.vector)
    ddply(df, "seqname", doByChr, covs, thresh)
}

if (!is.null(params$readsInput)) {
    # discard undercovered NFRs
    reads <- get(load(params$reads))
    nfr <- rmUnderCovered(nfr, reads, 0.75)
}

## Store output ###############################################################

message("-- saving gff output")
# save the output as a gff3 too
writeGff(df2gff(nfr,
                source="nucleR",
                feature="Nucleosme Free Region"),
         params[["output"]])

###############################################################################
