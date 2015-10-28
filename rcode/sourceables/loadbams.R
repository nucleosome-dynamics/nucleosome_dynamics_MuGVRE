#!/usr/bin/Rscript

# Functions to load reads from a BAM file.
# Can process either single-end or paired-end experiments

library(IRanges)
library(GenomicRanges)
library(Rsamtools)

SOURCE.DIR <- "/home/rilla/nucleServ/rcode/sourceables"
source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

loadSingleBam <- function(exp)
{
    # Rsamtools
    bam <- scanBam(exp,
                   param=ScanBamParam(what=c("pos",
                                             "qwidth",
                                             "strand",
                                             "rname")))[[1]]

    bam <- lapply(bam,
                  `[`,
                  Reduce(`&`,
                         lapply(bam,
                                function(x) !is.na(x))))

    # IRanges
    RangedData(space  = bam$rname,
               ranges = IRanges(start = bam[["pos"]],
                                width = bam[["qwidth"]]),
               strand = bam[["strand"]])
}

#Binary conversion
int2base <- function(x, b=2) {
    xi <- as.integer(x)
    if (any(is.na(xi) | ((x-xi)!=0))) {
        print(list(ERROR="x not integer", x=x))
    }
    N <- length(x)
    xMax <- max(x)
    ndigits <- (floor(logb(xMax, base=2)) + 1)
    Base.b <- array(NA, dim=c(N, ndigits))
    for (i in 1:ndigits) {
            Base.b[, ndigits-i+1] <- (x %% b)
            x <- (x %/% b)
    }
    if (N ==1) {
        Base.b[1, ]
    } else {
        Base.b
    }
}

bamFlagMatrix <- function(flags)
{
    bin <- int2base(flags)
    n <- ncol(bin)
    colnames(bin) <- c(rev(names(formals(scanBamFlag))[1:n]))
    return(bin)
}

processStrand <- function(strand, bam, flags)
{
    message(sprintf("processing strand %s", strand))
    is.paired <- flags[, "isPaired"] & flags[, "isProperPair"]
    mate1 <- flags[, "isFirstMateRead"]
    mate2 <- flags[, "isSecondMateRead"]
    strand.check <- flags[, "isMinusStrand"]

    if (strand == "+") {
        p1mate <- mate1
        p2mate <- mate2
    } else if (strand == "-") {
        p1mate <- mate2
        p2mate <- mate1
    }

    p1 <- vectorizedAll(is.paired, p1mate, !strand.check)
    p2 <- vectorizedAll(is.paired, p2mate, strand.check)

    reads1 <- lapply(bam, `[`, p1)
    reads2 <- lapply(bam, `[`, p2)

    #Consistency check
    test <- all(vectorizedAll(reads1$mpos  == reads2$pos,
                              reads2$mpos  == reads1$pos,
                              reads1$rname == reads2$rname))
    if (!test) {
        stop(sprintf("ERROR: Mate selection for %s strand is invalid",
                     strand))
    } else {
        RangedData(space  = as.character(reads1$rname),
                   ranges = IRanges(start = reads1$pos,
                                    end   = with(reads2,
                                                 pos + qwidth - 1)))
    }
}

loadPairedBam <- function(file)
{
    # Read BAM file (only one access to disk, intended for Shared Memory)
    message(sprintf("reading file %s", file))
    bam <- scanBam(file=file,
                   param=ScanBamParam(what=c("qname",
                                             "flag",
                                             "rname",
                                             "strand",
                                             "pos",
                                             "qwidth",
                                             "mrnm",
                                             "mpos")))[[1]]
    # We will process the flags in R
    # (an alternative is multiple scanBam calls...)
    message("processing flags")
    flags <- bamFlagMatrix(bam$flag)

    do.call(rbind,
            lapply(c("+", "-"),
                   processStrand,
                   bam,
                   flags))
}

loadBAM <- function(f, type="single")
{
    if (type == "single") {
        loadSingleBam(f)
    } else if (type == "paired") {
        loadPairedBam(f)
    } else {
        stop("type must be `single` or `paired`")
    }
}
