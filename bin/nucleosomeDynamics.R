#!/usr/bin/Rscript

# Position: region where the movement happens
# Type: change in the nucleosome map
# Score: magnitude of the change

# class: type of hotspot (see help for all possible types)
# nuc: to which nucleosome the movement belongs. NA means that the hostpot couldn't be unequivocally associated to one nucleosome.
# number_of_reads: number of reads involved in this movement
# hreads: number of reads involved in the movement relative to the number of reads present in the area. This value ranges from 0 to 1 and the closest it is to 1, the more significant the movement.

## Imports ####################################################################

library(getopt)
library(NucDyn)

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
sourced <- c("helperfuns", "gff_funs")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

### Parameters and Arguments ###################################################

spec <- matrix(c("input1",         "a", 1, "character",
                 "input2",         "b", 1, "character",
                 "outputGff",      "c", 1, "character",
                 "plotRData",      "d", 1, "character",
                 "dynRData",       "e", 1, "character",
                 "threshRData",    "f", 1, "character",
                 "cores",          "g", 1, "integer",
                 "maxLen",         "h", 1, "integer",
                 "equalSize",      "i", 1, "logical",
                 "roundPow",       "j", 1, "integer",
                 "readSize",       "k", 1, "integer",
                 "maxDiff",        "l", 1, "integer",
                 "rangeStart",     "m", 1, "integer",
                 "rangeEnd",       "n", 1, "integer",
                 "chr",            "o", 1, "chracter",
                 "nuc.width",      "p", 1, "integer",
                 "combined",       "q", 1, "logical",
                 "same.magnitude", "r", 1, "integer",
                 "threshold",      "s", 1, "character",
                 "rep1",           "t", 1, "character",
                 "rep2",           "u", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

defaults <- list(cores          = 1,
                 maxLen         = 170,
                 equalSize      = FALSE,
                 plotRData      = NULL,
                 dynRData       = NULL,
                 threshRData    = NULL,
                 roundPow       = 5,
                 readSize       = 140,
                 maxDiff        = NULL,
                 rangeStart     = NULL,
                 rangeEnd       = NULL,
                 chr            = NULL,
                 scale          = 2,
                 nuc.width      = 120,
                 combined       = TRUE,
                 same.magnitude = 2,
                 threshold      = "60%",
                 rep1           = NULL,
                 rep2           = NULL)

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

if (is.null(params$rangeStart) || is.null(params$rangeEnd)) {
    params$range <- c()
} else {
    params$range <- c(rangeStart, rangeEnd)
}

if (is.null(params$maxDiff)) {
    params$maxDiff <- params$readSize/2
}

## Pipeline Itself ############################################################

if (!is.null(params$dynRData) && file.exists(params$dynRData)) {
    hs <- get(load(params$dynRData))

} else {
    r1 <- get(load(params$input1))
    r2 <- get(load(params$input2))

    message("running NucleosomeDynamics")
    dyn <- nucleosomeDynamics(setA      = r1,
                              setB      = r2,
                              maxLen    = params$maxLen,
                              equalSize = params$equalSize,
                              roundPow  = params$roundPow,
                              readSize  = params$readSize,
                              maxDiff   = params$maxDiff,
                              mc.cores  = params$cores)

    if (!is.null(params$plotRData)) {
        plotable <- makePlotable(dyn)
        save(plotable, file=params$plotRData)
    }

    message("finding hotspots")
    hs <- findHotspots(dyn            = dyn,
                       range          = params$range,
                       chr            = params$chr,
                       nuc.width      = params$nuc.width,
                       combined       = FALSE,
                       same.magnitude = params$same.magnitude,
                       threshold      = NULL,
                       mc.cores       = params$cores)

    if (!is.null(params$dynRData)) {
        save(hs, file=params$dynRData)
    }
}

if (!is.null(params$rep1) && !is.null(params$rep1)) {
    if (grepl("%$", params$threshold)) {
        params$threshold <- defaults$scale
    }

    if (!is.null(params$threshRData) && file.exists(params$threshRData)) {
        thresh <- get(load(params$threshRData))
    } else {
        rep.fs <- strsplit(c(params$rep1, params$rep2), ",")
        pairs <- do.call(mapply, c(list, rep.fs, SIMPLIFY=FALSE))
        reps <- lapply(pairs, lapply, compose(get, load))
        thresh <- do.call(getVariableThreshold,
                          c(reps, mc.cores=params$cores))

        if (!is.null(params$threshRData)) {
            save(thresh, file=params$threshRData)
        }
    }
    thresh@scale <- as.numeric(params$threshold)
} else if (grepl("%$", params$threshold)) {
    thresh <- params$threshold
} else {
    thresh <- as.numeric(params$threshold)
}

hs <- applyThreshold(hs, thresh)

if (params$combined) {
    message("combining")
    hs <- combiner(hs,
                   params$nuc.width,
                   params$same.magnitude,
                   mc.cores=params$cores)
}

### Store the Result ###########################################################

#Score
#coord
#class
#nuc
#number_of_reads
#readsInvolved
#hreads

hs$nuc[hs$nuc == 0] <- NA
names(hs)[names(hs) == "type"] <- "class"
names(hs)[names(hs) == "chr"] <- "seqname"
names(hs)[names(hs) == "nreads"] <- "score"
#names(hs)[names(hs) == "totalReads"] <- "number_of_reads"
names(hs)[names(hs) == "readsInvolved"] <- "number_of_reads"

hs$coord <- NULL
hs$totalReads <- NULL
hs$freads <- NULL

message("saving output as gff")
writeGff(df2gff(hs,
                source="NucleosomeDynamics",
                feature="Nucleosome change"),
         params$outputGff)
message("done")

###############################################################################
