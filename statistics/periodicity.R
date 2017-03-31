#!/usr/bin/Rscript


## Imports ####################################################################

library(getopt)
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

spec <- matrix(c("input",     "a", 1, "character",
                 "genome",    "b", 1, "character"
                 "out_genes", "c", 1, "character",
                 "out_gw",    "d", 1, "character"),
               byrow=TRUE,
               ncol=4)

params <- getopt(spec)

## Read genes #################################################################

message("-- loading used genome")
genes <- getGenes(params$genome)

genes$tss <- as.numeric(genes$tss)
genes$tts <- as.numeric(genes$tts)

## Statistics per gene ########################################################

message("-- computing statistics per gene")
phase <- readGff(params$input)

genes_out = phase[, c("id", "score_phase", "score_autocorrelation")]
genes_out = merge(genes[, c("name", "tss")],
                  genes_out,
                  by.x="name",
                  by.y="id",
                  all.x=TRUE)
genes_out = genes_out[, c(1,3,4)]
genes_out = rbind(c("Name", "Score phase", "Score autocorrelation"), 
                  genes_out)

write.table(genes_out,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

gw_out = data.frame(c("Phased genes", "Not-phased genes", "Other genes"),
                    c(sum(phase$score_phase <= 25),
                      sum(phase$score_phase >= 56),
                      sum(phase$score_phase > 25 & phase$score_phase <= 55)))

write.table(gw_out,
            params$out_gw,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")
