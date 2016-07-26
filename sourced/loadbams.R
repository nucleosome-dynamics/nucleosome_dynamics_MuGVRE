#!/usr/bin/Rscript

# Functions to load reads from a BAM file.
# Can process either single-end or paired-end experiments

library(IRanges)
library(GenomicRanges)
library(Rsamtools)

source(paste(SOURCE.DIR,
             "helperfuns.R",
             sep="/"))

sortBy <- function (xs, a)
    lapply(xs, `[`, sort.list(xs[[a]]))

loadSingleBam <- function(exp)
{
    what <- c("pos", "qwidth", "strand", "rname")
    bam <- scanBam(exp, param=ScanBamParam(what=what))[[1]]

    non.na <- Reduce(`&`, lapply(bam, Negate(is.na)))
    filtered.bam <- lapply(bam, `[`, non.na)

    # IRanges
    RangedData(space  = filtered.bam$rname,
               ranges = IRanges(start = filtered.bam[["pos"]],
                                width = filtered.bam[["qwidth"]]),
               strand = filtered.bam[["strand"]])
}

#Binary conversion
int2base <- function(x, b=2)
{
    xi <- as.integer(x)
    if (any(is.na(xi) | ((x-xi) != 0))) {
        print(list(ERROR="x not integer", x=x))
    }
    N <- length(x)
    xMax <- max(x)
    ndigits <- floor(logb(xMax, base=2)) + 1
    Base.b <- array(NA, dim=c(N, ndigits))
    for (i in 1:ndigits) {
        Base.b[, ndigits-i+1] <- (x %% b)
        x <- (x %/% b)
    }
    if (N == 1) {
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

    #is.paired <- flags[, "isPaired"] & flags[, "isProperPair"]
    #mate1 <- flags[, "isFirstMateRead"]
    #mate2 <- flags[, "isSecondMateRead"]
    #strand.check <- flags[, "isMinusStrand"]

    #if (strand == "+") {
    #    p1mate <- mate1
    #    p2mate <- mate2
    #} else if (strand == "-") {
    #    p1mate <- mate2
    #    p2mate <- mate1
    #}

    ## Separate the paired reads
    #p1 <- vectorizedAll(is.paired, p1mate, !strand.check)
    #p2 <- vectorizedAll(is.paired, p2mate, strand.check)
    ##unsorted.reads1 <- lapply(bam, `[`, p1)
    ##unsorted.reads2 <- lapply(bam, `[`, p2)

    p1 <- ifelse(strand == "+", 99, 163)
    p2 <- ifelse(strand == "+", 147, 83)

    unsorted.reads1 <- bam[flags == p1, ]
    unsorted.reads2 <- bam[flags == p2, ]

    rownames(unsorted.reads1) <- as.vector(unsorted.reads1$qname)
    rownames(unsorted.reads2) <- as.vector(unsorted.reads2$qname)

    # Sort by the name of the reads. Assiming the paired reads will have the
    # same name, this will keep the pairs in the same position
    common <- intersect(rownames(unsorted.reads1), rownames(unsorted.reads2))

    reads1 <- unsorted.reads1[common, ]
    reads2 <- unsorted.reads2[common, ]

    # Consistency check
    test <- all(vectorizedAll(reads1$mpos  == reads2$pos,
                              reads2$mpos  == reads1$pos,
                              reads1$rname == reads2$rname))

    if (!test) {
        stop(sprintf("ERROR: Mate selection for %s strand is invalid",
                     strand))
    } else {
        RangedData(space  = as.character(reads1$rname),
                   ranges = IRanges(start = reads1$pos,
                                    end   = reads2$pos + reads2$qwidth - 1))
    }
}

loadPairedBam <- function(file)
{
    # Read BAM file (only one access to disk, intended for Shared Memory)
    message(sprintf("reading file %s", file))

    what <- c("qname",
              "flag",
              "rname",
              "strand",
              "pos",
              "qwidth",
              "mrnm",
              "mpos")
    bam <- as.data.frame(scanBam(file=file, param=ScanBamParam(what=what))[[1]])

    message("processing flags")
    ## We will process the flags in R
    ## (an alternative is multiple scanBam calls...)
    #flags <- bamFlagMatrix(bam$flag)
    flags <- bam$flag %% 256

    # Process both strand and return the reads in sorted order
    sortReads(rbind(processStrand("+", bam, flags),
                    processStrand("-", bam, flags)))
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
