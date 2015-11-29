#!/usr/bin/env Rscript

## Imports ####################################################################

library(IRanges)
library(parallel)
library(getopt)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
sourced <- c("helperfuns", "wig_funs", "get_genes", "periodicity_funs",
             "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

defaults <- list(period = 160)

spec <- matrix(c("calls",       "a", 1, "character",
                 "genome",      "c", 1, "character",
                 "output",      "d", 1, "character",
                 "cores",       "e", 1, "integer",
                 "periodicity", "f", 1, "double"),
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
genes <- readGenome(params$genome,
                    cols=c("TXSTART",
                           "TXEND",
                           "TXCHROM",
                           "TXSTRAND"))
message("reading calls")
calls.df <- readGff(params$calls)
calls.rd <- with(calls.df,
                 RangedData(space=seqname,
                            range=IRanges(start=start,
                                          end=end)))

## Do it ######################################################################

message("identifying first and last nucleosomes")
genes.nucs <- findGenesNucs(genes, calls.rd, params$mc.cores)
covPredAll <- getPeriodCov(genes.nucs, params$period, params$mc.cores)

## Store output ###############################################################

message("writting bigWig output")
splited <- lapply(covPredAll, splitAtZeros)
writeBigWig(splited, params$output, params$genome)

##############################################################################
