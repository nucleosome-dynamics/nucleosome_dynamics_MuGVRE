#!/usr/bin/Rscript

## Imports ####################################################################

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

###############################################################################

getFirstTx <- function (x, df)
{
    entries <- subset(df, GENEID == x)
    f <- `[[`(list("+"=which.min,
                   "-"=which.max),
              unique(entries$TXSTRAND))
    entries[f(entries$TXSTART), ]
}

cleanExons <- function (df)
{
    dupls <- myFilter(df$GENEID,
                      duplicated)
    sortDfBy(rbind(subset(df,
                          !GENEID %in% dupls),
                   do.call(rbind,
                           lapply(dupls,
                                  getFirstTx,
                                  df))),
             c("TXCHROM",
               "TXSTART"))
}

readGenome <- function (genome, cols)
{
    genome.lib <- paste0("TxDb.Scerevisiae.UCSC.",
                         genome,
                         ".sgdGene")
    library(genome.lib,
            character.only=TRUE)
    genome <- get(genome.lib)
    cleanExons(suppressWarnings(select(genome,
                                       keys=keys(genome),
                                       columns=cols,
                                       keytype="GENEID")))
}

###############################################################################
