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

spec <- matrix(c("input",     "a", 1, "character",
                 "genome",    "b", 1, "character",
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

genes_gr <- GRanges(genes$chrom,
                    IRanges(start=genes$start, end=genes$end),
                    name=genes$name)

## Statistics per gene ########################################################

message("-- computing statistics per gene")
nd <- readGff(params$input)

nd_gr <- GRanges(nd$seqname,
                 IRanges(start=nd$start, end=nd$end),
                 class=nd$class)

incl = nd_gr[nd_gr$class == "INCLUSION"]
evic = nd_gr[nd_gr$class == "EVICTION"]
shift_p = nd_gr[nd_gr$class == "SHIFT +"]
shift_m = nd_gr[nd_gr$class == "SHIFT -"]
incr_fuzzy = nd_gr[nd_gr$class == "INCREASED FUZZYNESS"]
decr_fuzzy = nd_gr[nd_gr$class == "DECREASED FUZZYNESS"]

genes$nIncl = countOverlaps(genes_gr, incl)
genes$nEvic = countOverlaps(genes_gr, evic)
genes$nShift_p = countOverlaps(genes_gr, shift_p)
genes$nShift_m = countOverlaps(genes_gr, shift_m)
genes$nIncr_fuzzy = countOverlaps(genes_gr, incr_fuzzy)
genes$nDecr_fuzzy = countOverlaps(genes_gr, decr_fuzzy)

i <- c("name",
       "nIncl",
       "nEvic",
       "nShift_p",
       "nShift_m",
       "nIncr_fuzzy",
       "nDecr_fuzzy")
stat_nd <- genes[, i]

stat_nd <- rbind(c("Name",
                   "Inclusions",
                   "Evictions",
                   "Shifts+",
                   "Shifts-",
                   "Increased Fuzziness",
                   "Decreased Fuzziness"),
                 stat_nd)

write.table(stat_nd,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

#--- Mean and std.dev ---
message("-- computing statistics genome-wide")

nd_tab = table(nd$class)/sum(table(nd$class))

i <- c("INCLUSION",
       "EVICTION",
       "SHIFT +",
       "SHIFT -",
       "INCREASED FUZZYNESS",
       "DECREASED FUZZYNESS")
nd_tab = nd_tab[i]
names(nd_tab) = c("Inclusion",
                  "Eviction",
                  "Shift +",
                  "Shift -",
                  "Increased Fuzziness",
                  "Decreased Fuzziness")

png(params$out_gw)
par(mar=c(6,4,2,2) + 0.1)
bp = barplot(nd_tab,
             col=c("#04B431","#FF0000","#8000FF","#0040FF","#424242","#BDBDBD"),
             xlab="",
             ylab="Proportion",
             xaxt="n")
text(x=bp, y=-0.05, names(nd_tab), xpd=TRUE, srt=45)
dev.off()
