#!/usr/bin/Rscript

library(NucDyn)
library(IRanges)
library(GenomicRanges)

in.dir <- "/orozco/services/R-data/tmp_wd"
in.f <- "120502_SN365_B_L002_GGM-34_120502_SN365_B_L002_GGM-35.RData"

dyn <- get(load(paste(in.dir, in.f, sep="/")))

params <- list(maxLen         = 170,
               equalSize      = FALSE,
               roundPow       = 5,
               readSize       = 140,
               maxDiff        = 70,
               rangeStart     = NULL,
               rangeEnd       = NULL,
               chr            = NULL,
               combined       = TRUE,
               same.magnitude = 2,
               threshold      = "60%",
               ARGS           = character(0),
               input1         = paste0("/orozco/services/R-data/tmp_wd",
                                       "/120502_SN365_B_L002_GGM-34.RData"),
               input2         = paste0("/orozco/services/R-data/tmp_wd",
                                       "/120502_SN365_B_L002_GGM-35.RData"),
               outputGff      = paste0("/orozco/services/R-data/tmp_wd",
                                       "/120502_SN365_B_L002_GGM-34",
                                       "_120502_SN365_B_L002_GGM-35.gff"),
               outputRData    = paste0("/orozco/services/R-data/tmp_wd",
                                       "/120502_SN365_B_L002_GGM-34_120502",
                                       "_SN365_B_L002_GGM-35.RData"),
               mc.cores       = 1)

message("finding hotspots")
hs <- findHotspots(dyn            = dyn,
                   #range          = params$range,
                   chr            = params$chr,
                   nuc.width      = params$nuc.width,
                   combined       = params$combined,
                   same.magnitude = params$same.magnitude,
                   threshold      = params$threshold,
                   mc.cores       = params$mc.cores)
