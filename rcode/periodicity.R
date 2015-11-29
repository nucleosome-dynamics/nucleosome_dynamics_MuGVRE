#!/usr/bin/env Rscript

library(IRanges)
library(parallel)
library(getopt)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

spec <- matrix(c("calls",             "a", 1, "character",
                 "genome",            "c", 1, "character",
                 "output",            "d", 1, "character",
                 "cores",             "e", 1, "integer"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

#genome.name <- "TxDb.Scerevisiae.UCSC.sacCer3.sgdGene"
#calls.f <- "/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.gff"
#out.f <- "/home/rilla/period.bw"

# Periodic positioning of nucleosomes from p1.pos and last.pos

period <- 165

# Create coverage given a vector and a period
ecov <- function(x, period)
    (1 + sin(pi/2 + 2*pi/period*x)) * 0.8^(abs(x)/period)

###############################################################################

cleanExons <- function (df)
{
    dupls <- myFilter(df$GENEID, duplicated)
    sortDfBy(rbind(subset(df,
                          !GENEID %in% dupls),
                   do.call(rbind,
                           lapply(dupls,
                                  getFirstTx,
                                  df))),
             c("TXCHROM", "TXSTART"))
}

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

splitAtZeros <- function (xs)
{
    xs <- as.vector(xs)
    counts <- rle(xs != 0)

    jdxs <- cumsum(counts$length)
    idxs <- c(1, jdxs[-length(jdxs)] + 1)

    by.kinds <- mapply(function(i, j) xs[i:j],
                       idxs, jdxs,
                       SIMPLIFY=FALSE)
    names(by.kinds) <- idxs
    return(by.kinds[counts$values])
}

writeWig <- function (x, outf)
{
    tag.rows <- unlist(lapply(names(x),
                              function(chr)
                                  paste("fixedStep",
                                        paste0("chrom=", chr),
                                        paste0("start=", names(x[[chr]])),
                                        "step=1",
                                        sep=" ")))
    vals <- unlist(x, recursive=FALSE)
    idx <- order(c(seq_along(tag.rows), seq_along(vals)))
    cat(unlist(c(tag.rows, vals)[idx]),
        sep="\n",
        file=outf)
}

message("loading genome")
library(args$genome, character.only=TRUE)
genome <- get(args$genom)

genes <- cleanExons(suppressWarnings(select(genome,
                                            keys=keys(genome),
                                            columns=c("TXSTART",
                                                      "TXEND",
                                                      "TXCHROM",
                                                      "TXSTRAND"),
                                            keytype="GENEID")))

message("reading calls")
calls.df <- readGff(args$calls)
calls.rd <- with(calls.df,
                 RangedData(space=seqname,
                            range=IRanges(start=start,
                                          end=end)))

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
covPredAll <- mclapply(
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
                  period)),
    mc.cores=1
)
names(covPredAll) <- chroms

message("writting wig output")
splited <- lapply(covPredAll, splitAtZeros)
wigf <- sub(".bw$", ".wig", args$output)
writeWig(splited, wigf)
