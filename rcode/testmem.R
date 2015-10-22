#!/opt/R-2.15.3/bin/Rscript

## Imports ####################################################################

library(getopt)
library(NucDyn)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))
source(paste(SOURCE.DIR,
             "loadbams.R",
             sep="/"))

## Parameters and Arguments ###################################################

defaults <- list(maxLen         = 170,
                 equalSize      = FALSE,
                 roundPow       = 5,
                 readSize       = 140,
                 maxDiff        = NULL,
                 rangeStart     = NULL,
                 rangeEnd       = NULL,
                 chr            = NULL,
                 combined       = TRUE,
                 same.magnitude = 2,
                 threshold      = "60%")

spec <- matrix(c("input1",         "a", 1, "character",
                 "input2",         "b", 1, "character",
                 "output",         "c", 1, "character",
                 "cores",          "d", 1, "integer",
                 "maxLen",         "e", 1, "integer",
                 "equalSize",      "f", 1, "logical",
                 "roundPow",       "g", 1, "integer",
                 "readSize",       "h", 1, "integer",
                 "maxDiff",        "i", 1, "integer",
                 "rangeStart",     "j", 1, "integer",
                 "rangeEnd",       "k", 1, "integer",
                 "chr",            "l", 1, "chracter",
                 "combined",       "m", 1, "logical",
                 "same.magnitude", "n", 1, "integer",
                 "threshold",      "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- sub("cores", "mc.cores", names(args))

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

#> #> #> #> #> #>


#> #> #> #> #> #>

## Pipeline Itself ############################################################

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
                          mc.cores  = params$mc.cores)

message("finding hotspots")
hs <- findHotspots(dyn            = dyn,
                   range          = params$range,
                   chr            = params$chr,
                   nuc.width      = params$nuc.width,
                   combined       = params$combined,
                   same.magnitude = params$same.magnitude,
                   threshold      = params$threshold,
                   mc.cores       = params$mc.cores)

## Store the Result ###########################################################

names(hs) <- sub("chr", "seqname", names(hs))
names(hs) <- sub("nreads", "score", names(hs))

message("saving output as gff")
writeGff(df2gff(hs,
                source="NucleosomeDynamics",
                feature="Nucleosome change"),
         params$output)
message("done")

###############################################################################
