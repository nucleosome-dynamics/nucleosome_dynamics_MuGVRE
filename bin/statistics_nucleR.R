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

genes_gr <- GRanges(genes$chrom, IRanges(start=genes$start, end=genes$end), name=genes$name)


## Statistics per gene ########################################################

message("-- computing statistics per gene")
nuc <- readGff(params$input)


nuc_gr <- GRanges(nuc$seqname, IRanges(start=nuc$start, end=nuc$end), class=nuc$class)

nuc_well <- nuc_gr[nuc_gr$class=="W"]
nuc_fuzzy <- nuc_gr[nuc_gr$class=="F"]
nuc_uncertain <- nuc_gr[nuc_gr$class=="uncertain"]

genes$TotalNucleosomes <- countOverlaps(genes_gr, nuc_gr)
genes$TotalWellPositioned <- countOverlaps(genes_gr, nuc_well)
genes$TotalFuzzy <- countOverlaps(genes_gr, nuc_fuzzy)
genes$TotalUncertain <- countOverlaps(genes_gr, nuc_uncertain)

genes_out = genes[,c("name", "TotalNucleosomes", "TotalWellPositioned", "TotalFuzzy", "TotalUncertain")]
genes_out = rbind(c("Name", "Total Nucleosomes", "Total Well-Positioned", "Total Fuzzy", "Total Uncertain"), 
              genes_out)


OUT_GENES =  gsub(".gff", "_genes_stats.csv", params$input)
tmp = strsplit(OUT_GENES, "/")[[1]][length(strsplit(OUT_GENES, "/")[[1]])]

OUT_GENES = gsub(tmp, paste(".", tmp, sep=""), OUT_GENES, fixed=T)

write.table(genes_out, OUT_GENES, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")


## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")
nuc_gr$class = as.character(nuc_gr$class)
class_lab = gsub("F", "Fuzzy", nuc_gr$class)
class_lab = gsub("W", "Well-positioned", class_lab)
class_lab = gsub("uncertain", "Uncertain", class_lab)

gw_stat <- as.data.frame(table(class_lab))
colnames(gw_stat) = c("Class", "Frequency")
gw_stat$Class = as.character(gw_stat$Class)
gw_stat <- rbind(gw_stat, c("Total", length(nuc_gr)))


OUT_GW <- gsub(".gff", "_stats.csv", params$input)
tmp = strsplit(OUT_GW, "/")[[1]][length(strsplit(OUT_GW, "/")[[1]])]

OUT_GW = gsub(tmp, paste(".", tmp, sep=""), OUT_GW, fixed=T)

write.csv(gw_stat, OUT_GW, row.names=F, quote=F)



