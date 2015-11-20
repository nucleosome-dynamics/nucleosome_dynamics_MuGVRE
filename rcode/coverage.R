#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(htSeqTools)
library(nucleR)
library(IRanges)
library(GenomicRanges)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

tobig.bin <- "/home/rilla/nucleServ/wig_utils/wigToBigWig"
chrom.sizes.f <- "/home/rilla/nucleServ/wig_utils/sacCer3.chrom.sizes"

#params$output <- "/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.bw"
#params$mc.cores <- 1
#params$input <- "/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
#params$type <- "paired"


## Parameters and Arguments ###################################################

defaults <- list(type           = NULL,
                 fdrOverAmp     = 0.05,
                 components     = 1,
                 fragmentLen    = NULL,
                 trim           = 50,

spec <- matrix(c("input",       "a", 1, "character",
                 "output",      "b", 1, "character",
                 "cores",       "c", 1, "integer",
                 "type",        "e", 1, "character",
                 "fdrOverAmp",  "d", 1, "double",
                 "components",  "f", 1, "integer",
                 "fragmentLen", "g", 1, "integer",
                 "trim",        "h", 1, "integer")
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

if (!grepl(".bw$", params$output)) {
    message("Expected output is a big wig file")
    q("no")
}

## Pipeline Itself ############################################################

reads <- get(load(params$input))

message("filtering duplicated reads")
f.reads <- filterDuplReads(reads,
                           fdrOverAmp=params$fdrOverAmp,
                           components=params$components,
                           mc.cores=params$mc.cores)

if (is.null(params$fragmentLen)) {
    if (params$type == "single") {
        message("estimating fragment length")
        params$fragmentLen <- fragmentLenDetect(f.reads)
    } else if (params$type == "paired") {
        params$fragmentLen <- 170
    }
}

message("processing reads")
prep <- processReads(f.reads,
                     type=params$type,
                     fragmentLen=params$fragmentLen,
                     trim=params$trim)

message("calculating coverage")
cover <- coverage.rpm(prep)

## Store the Result ###########################################################

# we'll only save non-zero regions in our wig
splited <- lapply(cover, splitAtZeros)

wigf <- sub(".bw$", ".wig", params$output)

message("writing wig file")
writeWig(splited, wigf)
message("converting to bigWig")
system(paste(tobig.bin, wigf, chrom.sizes.f, params$out))
file.remove(wigf)

###############################################################################
