#!/opt/R-2.15.3/bin/Rscript

library(getopt)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "loadbams.R",
             sep="/"))

spec <- matrix(c("type",   "t", 1, "character",
                 "input",  "i", 1, "character",
                 "output", "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

message("-- loading ", args[["input"]])
reads <- loadBAM(args[["input"]], args[["type"]])

gc()

message("-- saving ", args[["output"]])
save(reads, file=args[["output"]])
