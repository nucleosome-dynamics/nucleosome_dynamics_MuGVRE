#!/usr/bin/Rscript


## Imports ####################################################################

library(getopt)
library(htSeqTools)
library(nucleR)
library(IRanges)
library(GenomicRanges)

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
sourced <- c("get_genes", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}


## Parameters and Arguments ###################################################

defaults <- list()

spec <- matrix(c("input", "a", 1, "character",
                 "genome", "b", 1, "character"
                ),
               byrow=TRUE,
               ncol=4)

args <- getopt(spec)

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}




## Read NFR data  ########################################################


nfr <- readGff(params$input)
nfr_width = nfr$end - nfr$start


## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

stat_nfr = data.frame(NFR=c("Total", "Mean width", "Std. Dev. width"), 
                      Value= c(nrow(nfr), round(mean(nfr_width),2), round(sd(nfr_width),2)))


OUT_GW <- gsub(".gff", "_stats.csv", params$input)
tmp = strsplit(OUT_GW, "/")[[1]][length(strsplit(OUT_GW, "/")[[1]])]

OUT_GW = gsub(tmp, paste(".", tmp, sep=""), OUT_GW, fixed=T)

write.csv(stat_nfr, OUT_GW, row.names=F, quote=F)



