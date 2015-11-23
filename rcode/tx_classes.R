#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)
library(GenomicRanges)
library(htSeqTools)
library(nucleR)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))
source(paste(SOURCE.DIR,
             "nucleosome_patterns.R",
             sep="/"))

## Parameters and Arguments ###################################################

defaults <- list(window            = 300,
                 p1.max.merge      = 3,
                 p1.max.downstream = 20,
                 open.thresh       = 215,
                 max.uncovered     = 150)

spec <- matrix(c("calls",             "a", 1, "character",
                 "coverage",          "b", 1, "character",
                 "genome",            "c", 1, "character",
                 "output",            "d", 1, "character",
                 "cores",             "e", 1, "integer",
                 "window",            "f", 1, "integer",
                 "p1.max.merge",      "g", 1, "integer",
                 "p1.max.downstream", "h", 1, "integer",
                 "open.thresh",       "i", 1, "integer",
                 "max.uncovered",     "j", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some function declarations #################################################

getFirstTx <- function(x, df)
{
    entries <- subset(df, GENEID == x)
    f <- `[[`(list("+"=which.min,
                   "-"=which.max),
              unique(entries$TXSTRAND))
    entries[f(entries$TXSTART), ]
}

cleanExons <- function(df)
{
    dupls <- myFilter(df$GENEID,
                      duplicated)
    sortDfBy(rbind(subset(df,
                          !GENEID %in% dupls),
                   do.call(rbind,
                           lapply(dupls,
                                  getFirstTx,
                                  df))),
             c("TXCHROM",
               "TXSTART"))
}

## Load the inputs ############################################################

message("-- loading inputs")
nucs <- with(readGff(params$calls),
             RangedData(ranges  = IRanges(start = as.numeric(start),
                                          end   = as.numeric(end)),
                        space   = as.character(seqname),
                        score   = as.numeric(score),
                        score_w = as.numeric(score_w),
                        score_h = as.numeric(score_h),
                        nmerge  = as.numeric(nmerge),
                        class   = as.character(class)))

cover <- get(load(params$coverage))

message("-- loading used genome")
library(params$genome, character.only=TRUE)
genome <- get(params$genome)

genes <- cleanExons(select(genome,
                           keys=keys(genome),
                           columns=c("TXSTART", "TXCHROM", "TXSTRAND"),
                           keytype="GENEID"))

## Do it ######################################################################

message("-- checking the classes")
tx.classes <- with(params,
                   nucleosomePatternsDF(calls             = nucs,
                                        cover             = cover,
                                        df                = genes,
                                        col.id            = "GENEID",
                                        col.chrom         = "TXCHROM",
                                        col.pos           = "TXSTART",
                                        col.strand        = "TXSTRAND",
                                        window            = window,
                                        p1.max.merge      = p1.max.merge,
                                        p1.max.downstream = p1.max.downstream,
                                        open.thresh       = open.thresh,
                                        max.uncovered     = max.uncovered,
                                        mc.cores          = mc.cores))

## Store output ###############################################################

message("-- saving gff output")
names(tx.classes) <- subMany(c("chrom",   "pos"),
                             c("seqname", "start"),
                             names(tx.classes))
gff <- df2gff(tx.classes,
              source="nucleR",
              feature="txStart class",
              end=tx.classes$start)
writeGff(gff, params$output)

###############################################################################
