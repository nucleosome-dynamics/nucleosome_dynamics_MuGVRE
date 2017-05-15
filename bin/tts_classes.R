#!/usr/bin/env Rscript

library(getopt)
library(IRanges)
library(GenomicRanges)
library(htSeqTools)
library(nucleR)

where <- function() "/home/rilla/nucleServ/sourced"
SOURCE.DIR <- paste(where(), "../sourced", sep="/")
sourced <- c("helperfuns", "nucleosome_patterns", "get_genes", "gff_funs")
for (x in paste0(SOURCE.DIR, "/", sourced, ".R")) {
    source(x)
}

## Parameters and Arguments ###################################################

defaults <- list(mc.cores          = 1,
                 window            = 300,
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

nucs <- with(readGff(params$calls),
             RangedData(ranges  = IRanges(start = as.numeric(start),
                                          end   = as.numeric(end)),
                        space   = as.character(seqname),
                        score   = as.numeric(score),
                        score_w = as.numeric(score_w),
                        score_h = as.numeric(score_h),
                        nmerge  = as.numeric(nmerge),
                        class   = as.character(class)))

genes <- getGenes(params$genome)
genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

tss.classes <- with(params,
                    patternsByChrDF(calls             = nucs,
                                    df                = genes,
                                    col.id            = "name",
                                    col.chrom         = "chrom",
                                    col.pos           = "tss",
                                    col.strand        = "strand",
                                    window            = window,
                                    p1.max.merge      = p1.max.merge,
                                    p1.max.downstream = p1.max.downstream,
                                    open.thresh       = open.thresh,
                                    max.uncovered     = max.uncovered,
                                    position          = "tss",
                                    mc.cores          = 10))

tts.classes <- with(params,
                    patternsByChrDF(calls             = nucs,
                                    df                = genes,
                                    col.id            = "name",
                                    col.chrom         = "chrom",
                                    col.pos           = "tts",
                                    col.strand        = "strand",
                                    window            = window,
                                    p1.max.merge      = p1.max.merge,
                                    p1.max.downstream = p1.max.downstream,
                                    open.thresh       = open.thresh,
                                    max.uncovered     = max.uncovered,
                                    position          = "tts",
                                    mc.cores          = mc.cores))

gg <- subset(genes, chrom == "chrI")
tts <- subset(tts.classes, chrom == "chrI")
tss <- subset(tss.classes, chrom == "chrI")

rownames(tts) <- tts$id
rownames(tss) <- tss$id

tts <- tts[gg$name, ]
tss <- tss[gg$name, ]

gg <- gg[, c("name", "tss", "tts", "strand")]
tss <- tss[, c("id", "pos", "p1.pos", "m1.pos", "strand")]
tts <- tts[, c("id", "pos", "p1.pos", "m1.pos", "strand")]

gg[1, ]
tss[1, ]
tts[1, ]
