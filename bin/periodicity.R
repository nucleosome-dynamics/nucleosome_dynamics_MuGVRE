#!/usr/bin/env Rscript

## Imports ####################################################################

library(IRanges)
library(parallel)
library(getopt)
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
SOURCE.DIR <- "/home/rilla/nucleServ/sourced/"
sourced <- c("helperfuns", "wig_funs", "get_genes", "periodicity_funs",
             "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

#defaults <- list(period = 160)
defaults <- list(period   = 165,
                 mc.cores = 1,
                 genome = "R64-1-1")

spec <- matrix(c("calls",       "a", 1, "character",
                 "genome",      "c", 1, "character",
                 "bwOutput",    "d", 1, "character",
                 "gffOutput",   "h", 1, "character",
                 "cores",       "e", 1, "integer",
                 "coverage",    "f", 1, "character",
                 "periodicity", "g", 1, "double"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- subMany("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some function definitions ##################################################

message("loading genes")
genes <- getGenes(params$genome)
message("reading calls")
calls.df <- readGff(params$calls)
calls.rd <- with(calls.df,
                 RangedData(space=seqname,
                            range=IRanges(start=start,
                                          end=end)))

## Do it ######################################################################

message("identifying first and last nucleosomes")
cov <- get(load(params$coverage))

genes.nucs <- findGenesNucs(genes, calls.rd, params$mc.cores)
genes.nucs$dfi <- getDfi(genes.nucs$nuc.len, params$period)
genes.nucs$autocor <- autocorFromDf(genes.nucs, cov, params$period)

covPredAll <- getPeriodCov(genes.nucs, params$period, params$mc.cores)

## Store output ###############################################################

message("writting GFF")
names(genes.nucs)[names(genes.nucs) == "chrom"] <- "seqname"
names(genes.nucs)[names(genes.nucs) == "dfi"] <- "score"
gff <- df2gff(genes.nucs,
              source="nucleR",
              feature="nucleosome periodicity")
writeGff(gff, params$gffOutput)

message("writting bigWig output")
splited <- lapply(covPredAll, splitAtZeros)
writeBigWig(splited, params$bwOutput, params$genome)

##############################################################################
