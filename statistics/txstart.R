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
                 "genome",    "b", 1, "character",
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

## Statistics per gene ########################################################

message("-- computing statistics per gene")
tss <- readGff(params$input)

genes_out = tss[, c("id", "classification", "distance")]
genes_out = merge(genes[, c("name", "tss")],
                  genes_out,
                  by.x="name",
                  by.y="id",
                  all.x=TRUE)
genes_out = genes_out[, c(1, 3, 4)]
genes_out = rbind(c("Name", "TSS class", "Distance from -1 to +1"),
                  genes_out)

write.table(genes_out,
            params$out_genes,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE,
            sep=",")

## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

#--- Plot 1 ---

png(params$out_gw)
par(mar=c(c(5, 10, 4, 4) + 0.1))
tab_tss = table(tss$classification)
bpres = barplot(tab_tss, horiz=T,las=2, col="#66A61E")
text(x= tab_tss+80, y= bpres, labels=as.character(tab_tss), xpd=TRUE)
dev.off()


#--- Plot 2 ---

png(params$out_gw2)
plot(density(as.numeric(tss$distance), na.rm=T), 
     main="Distribution of NFR width around TSS", 
     xlab="Width",
     col="#1B9E77", lwd=2, ylim=c(0,0.01)
)
dev.off()
