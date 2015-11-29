#!/usr/bin/env Rscript

## Imports ####################################################################

library(IRanges)
library(parallel)
library(getopt)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
sourced <- c("helperfuns", "wig_funs", "get_genes")
for (x in sourced) {
    source(paste0(SOURCE.DIR, "/", x, ".R"))
}

## Parameters and Arguments ###################################################

defaults <- list(period = 160)

spec <- matrix(c("calls",       "a", 1, "character",
                 "genome",      "c", 1, "character",
                 "output",      "d", 1, "character",
                 "cores",       "e", 1, "integer",
                 "periodicity", "f", 1, "double"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

names(args) <- subMany("cores", "mc.cores", names(args))

params <- defaults
for (i in names(args)) {
    params[[i]] <- args[[i]]
}

## Some function definitions ##################################################

findFirstAndLast <- function(id, start, end, strand, dyads, chr)
{
    in.gene <- dyads[dyads >= start & dyads <= end]
    if (strand == "-") {
        in.gene <- rev(in.gene)
    }
    if (length(in.gene)) {
        data.frame(id=id,
                   chrom=chr,
                   strand=strand,
                   first=in.gene[1],
                   last=in.gene[length(in.gene)])
    } else {
        NULL
    }
}

ecov <- function (x, period)
    (1 + sin(pi/2 + 2*pi/period*x)) * 0.8^(abs(x)/period)

# Create a coverage for a gene given p1.pos and last.pos and their periodicity
coverageChr <- function(nuc.start, nuc.end, nuc.length, strand, L, period)
{
    revIt <- function(xs)
        rev(xs) * (-1)

    cov <- rep(0, L)

    for (i in seq_along(nuc.start)) {
        sper <- floor(period/2)
        xs <- (-sper):(period*floor(nuc.length[i]/period) + sper)

        if (strand[i] == "+") {
            a <- xs
            b <- revIt(xs)
        } else if (strand[i] == "-") {
            a <- revIt(xs)
            b <- xs
        }

        x <- nuc.start[i] + a
        y <- nuc.end[i] + b

        cov[round(x)] <- cov[round(x)] + ecov(a, period)
        cov[round(y)] <- cov[round(y)] + ecov(b, period)
    }

    cov
}

## Load inputs ################################################################

message("loading genes")
genes <- readGenome(params$genome,
                    cols=c("TXSTART",
                           "TXEND",
                           "TXCHROM",
                           "TXSTRAND"))
message("reading calls")
calls.df <- readGff(params$calls)
calls.rd <- with(calls.df,
                 RangedData(space=seqname,
                            range=IRanges(start=start,
                                          end=end)))

## Do it ######################################################################

chroms <- unique(genes$TXCHROM)
message("identifying first and last nucleosomes")
genes.by.chr <- lapply(chroms,
                       function (chr) subset(genes,
                                             TXCHROM == chr))
names(genes.by.chr) <- chroms
dyads.by.chr <- lapply(ranges(calls.rd), dyadPos)

chr.lens <- sapply(genes.by.chr,
                   compose(max,
                           partial(`[[`,
                                   "TXEND")))

used.cols <- c("GENEID", "TXSTART", "TXEND", "TXSTRAND")
genes.nucs <- do.call(rbind, lapply(
    chroms,
    function(chr)
        do.call(rbind,
                myFilter(iterDf(genes.by.chr[[chr]][, used.cols],
                                findFirstAndLast,
                                dyads.by.chr[[chr]],
                                chr),
                         Negate(is.null)))
))

genes.nucs$nuc.length <- abs(genes.nucs$last - genes.nucs$first)

# Compute Coverage for all genome
message("calculating theoretic coverage")
covPredAll <- xlapply(
    chroms,
    function (chr)
        do.call(coverageChr,
                c(unname(as.list(subset(genes.nucs,
                                        subset=chrom == chr,
                                        select=c("first",
                                                 "last",
                                                 "nuc.length",
                                                 "strand")))),
                  chr.lens[[chr]] + 500,
                  params$period)),
    mc.cores=params$mc.cores
)
names(covPredAll) <- chroms

## Store output ###############################################################

message("writting bigWig output")
splited <- lapply(covPredAll, splitAtZeros)
writeBigWig(splited, params$output, params$genome)

###############################################################################
