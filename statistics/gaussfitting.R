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
                 "out_gw",    "d", 1, "character",
                 "out_gw2",   "e", 1, "character"),
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
stf <- readGff(params$input)

stf_gr <- GRanges(stf$seqname,
                  IRanges(start=stf$start, end=stf$end),
                  score=stf$score)

ovlps <- findOverlaps(genes_gr, stf_gr)
split_ovlps <- split(subjectHits(ovlps), queryHits(ovlps))

res <- lapply(
    names(split_ovlps),
    function(i) {
        tmp = stf_gr[split_ovlps[[i]]]
        return(c(genes_gr$name[as.numeric(i)],
                 mean(tmp$score, na.rm=TRUE),
                 sd(tmp$score, na.rm=TRUE)))
    }
)

res <- do.call(rbind, res)
colnames(res) = c("name", "Mean_STF", "StdDev_STF")

stat_stf <- merge(genes, res, by="name", all.x=T)
colnames(stat_stf)[colnames(stat_stf)=="name"] = "Name"

write.table(stat_stf[, c("Name","Mean_STF","StdDev_STF")],
            params$out_genes,
            OUT_GENES,
            row.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

#--- Mean and std.dev ---
message("-- computing statistics genome-wide")

gw_stat = cbind(c("Mean stiffness", "Std. Dev. stiffness"),
                round(c(mean(stf_gr$score, na.rm=T),
                        sd(stf_gr$score, na.rm=T)),
                      4))

write.table(gw_stat,
            params$out_gw,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

#--- Plot distribution ---

stf$class=NA
stf$class[stf$score<0.1] = "0 - 0.1"
stf$class[stf$score>=0.1 & stf$score<0.2] = "0.1 - 0.2"
stf$class[stf$score>=0.2 & stf$score<0.3] = "0.2 - 0.3"
stf$class[stf$score>=0.3 & stf$score<0.4] = "0.3 - 0.4"
stf$class[stf$score>=0.4] = "0.4 - 1"

png(params$out_gw2)
barplot(table(stf$class)/sum(table(stf$class)),
        col=c("#98F5FF", "#71DAE2", "#4CC0C4", "#26A5A8", "#008B8B"),
        xlab="Stiffness", ylab="Proportion of genes")
dev.off()
