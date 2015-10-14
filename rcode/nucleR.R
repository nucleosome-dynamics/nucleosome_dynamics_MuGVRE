#!/opt/R-2.15.3/bin/Rscript

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
source(paste(SOURCE.DIR,
             "loadbams.R",
             sep="/"))

## Parameters and Arguments ###################################################

defaults <- list(type           = NULL,
                 fdrOverAmp     = 0.05,
                 components     = 1,
                 fragmentLen    = NULL,
                 trim           = 50,
                 pcKeepComp     = 0.02,
                 width          = 125,
                 threshold      = "35%",
                 dyad.length    = NULL,
                 min.overlap    = NULL,
                 score_w.thresh = 0.6,
                 score_h.thresh = 0.4)

spec <- matrix(c("input",       "a", 1, "character",
                 "output",      "b", 1, "character",
                 "cores",       "c", 1, "integer",
                 "type",        "e", 1, "character",
                 "fdrOverAmp",  "d", 1, "double",
                 "components",  "f", 1, "integer",
                 "fragmentLen", "g", 1, "integer",
                 "trim",        "h", 1, "integer",
                 "pcKeepComp",  "i", 1, "double",
                 "width",       "j", 1, "integer",
                 "threshold",   "k", 1, "character",
                 "dyadlength",  "l", 1, "integer",
                 "minoverlap",  "m", 1, "integer",
                 "wthresh",     "n", 1, "double",
                 "hthresh",     "o", 1, "double"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))
names(args) <- sub("dyadlength", "dyad.length", names(args))
names(args) <- sub("minoverlap", "min.overlap", names(args))
names(args) <- sub("wthresh", "score_w.thresh", names(args))
names(args) <- sub("hthresh", "score_h.thresh", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

if (is.null(params$dyad.length)) {
    params$dyad.length <- params$trim
}
if (is.null(params$min.overlap)) {
    params$min.overlap <- params$trim
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

message("filtering with filterFFT")
fft <- filterFFT(cover,
                 pcKeepComp=params$pcKeepComp,
                 mc.cores=params$mc.cores)

message("detecting peaks")
peaks <- peakDetection(fft,
                       width=params$width,
                       threshold=params$threshold,
                       score=FALSE,
                       mc.cores=params$mc.cores)

message("scoring peaks")
scores <- peakScoring(peaks,
                      fft,
                      threshold=params$threshold,
                      dyad.length=params$dyad.length,
                      mc.cores=params$mc.cores)

message("merging peaks")
merged <- mergeCalls(scores,
                     min.overlap=params$min.overlap,
                     mc.cores=params$mc.cores)

merged$type <- getType(merged$score_w,
                       merged$score_h,
                       params$score_w.thresh,
                       params$score_h.thresh)

## Store the Result ###########################################################

message("saving output as gff")
writeGff(df2gff(rd2df(merged),
                source="nucleR",
                feature="Nucleosome"),
         params$output)

message("done")

###############################################################################
