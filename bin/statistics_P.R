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


## Read genes #################################################################

message("-- loading used genome")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)


## Statistics per gene ########################################################

message("-- computing statistics per gene")
phase <- readGff(params$input)

genes_out = phase[,c("id", "score_phase", "score_autocorrelation")]
genes_out = merge(genes[, c("name", "tss")], genes_out, by.x="name", by.y="id", all.x=T)
genes_out = genes_out[,c(1,3,4)]
genes_out = rbind(c("Name", "Score phase", "Score autocorrelation"), 
                  genes_out)


OUT_GENES =  gsub(".gff", "_genes_stats.csv", params$input)
tmp = strsplit(OUT_GENES, "/")[[1]][length(strsplit(OUT_GENES, "/")[[1]])]

OUT_GENES = gsub(tmp, paste(".", tmp, sep=""), OUT_GENES, fixed=T)

write.table(genes_out, OUT_GENES, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")


## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

gw_out = data.frame(c("Phased genes", "Not-phased genes", "Other genes"), 
                    c(sum(phase$score_phase<=25), 
                      sum(phase$score_phase>=56),
                      sum(phase$score_phase>25 & phase$score_phase<=55)
                   ))

OUT_GW <- gsub(".gff", "_stats.csv", params$input)
tmp = strsplit(OUT_GW, "/")[[1]][length(strsplit(OUT_GW, "/")[[1]])]

OUT_GW = gsub(tmp, paste(".", tmp, sep=""), OUT_GW, fixed=T)


write.table(gw_out, OUT_GW, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")

