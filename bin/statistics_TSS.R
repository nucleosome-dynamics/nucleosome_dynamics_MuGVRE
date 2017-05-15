#!/usr/bin/Rscript


## Imports ####################################################################

library(getopt)
#library(htSeqTools)
#library(nucleR)
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
tss <- readGff(params$input)

genes_out = tss[,c("id", "classification", "distance")]
genes_out = merge(genes[, c("name", "tss")], genes_out, by.x="name", by.y="id", all.x=T)
genes_out = genes_out[,c(1,3,4)]
genes_out = rbind(c("Name", "TSS class", "Distance from -1 to +1"), 
              genes_out)


OUT_GENES =  gsub(".gff", "_genes_stats.csv", params$input)
tmp = strsplit(OUT_GENES, "/")[[1]][length(strsplit(OUT_GENES, "/")[[1]])]

OUT_GENES = gsub(tmp, paste(".", tmp, sep=""), OUT_GENES, fixed=T)

write.table(genes_out, OUT_GENES, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")


## Statistics genome-wide  ####################################################

message("-- computing statistics genome-wide")

#--- Plot 1 ---
OUT_GW <- gsub(".gff", "_stats1.png", params$input)
tmp = strsplit(OUT_GW, "/")[[1]][length(strsplit(OUT_GW, "/")[[1]])]

OUT_GW = gsub(tmp, paste(".", tmp, sep=""), OUT_GW, fixed=T)

png(OUT_GW)
par(mar=c(c(5, 10, 4, 4) + 0.1))
tab_tss = table(tss$classification)
bpres = barplot(tab_tss, horiz=T,las=2, col="#66A61E")
text(x= tab_tss+80, y= bpres, labels=as.character(tab_tss), xpd=TRUE)
dev.off()


#--- Plot 2 ---
OUT_GW2 <- gsub(".gff", "_stats2.png", params$input)
tmp = strsplit(OUT_GW2, "/")[[1]][length(strsplit(OUT_GW2, "/")[[1]])]

OUT_GW2 = gsub(tmp, paste(".", tmp, sep=""), OUT_GW2)


png(OUT_GW2)
plot(density(as.numeric(tss$distance), na.rm=T), 
     main="Distribution of NFR width around TSS", 
     xlab="Width",
     col="#1B9E77", lwd=2, ylim=c(0,0.01)
)
dev.off()



