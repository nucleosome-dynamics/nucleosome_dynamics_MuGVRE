#!/opt/R-2.15.3/bin/Rscript

library(getopt)
SOURCE.DIR <- paste0("/orozco/scratch/xesh0/orozco/rilla/nucleosome_dynamics",
                     "/webserver/sourced_funs")
source(paste(SOURCE.DIR, "loadbams.R", sep="/"))

spec <- matrix(c("type",   "t", 1, "character",
                 "input",  "i", 1, "character",
                 "output", "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

message("-- loading ", args[["input"]])
reads <- loadBAM(args[["input"]], args[["type"]])

message("-- saving ", args[["output"]])
save(reads, file=args[["output"]])
