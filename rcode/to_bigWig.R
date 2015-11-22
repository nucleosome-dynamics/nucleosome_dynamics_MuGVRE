#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)

tobig.bin <- "/home/rilla/nucleServ/wig_utils/wigToBigWig"

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

writeWig <- function(x, outf)
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

## Parameters and Arguments ###################################################

defaults <- list(type           = NULL,
                 fdrOverAmp     = 0.05,
                 components     = 1,
                 fragmentLen    = NULL,
                 trim           = 50)

spec <- matrix(c("input",  "i", 1, "character",
                 "output", "o", 1, "character",
                 "genome", "g", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

#args$input <- "/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.RData"
#args$output <- "/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.bw"
#args$genome <- "sacCer3"

###############################################################################

cover <- get(load(args$input))

# we'll only save non-zero regions in our wig
splited <- lapply(cover, splitAtZeros)

wigf <- sub(".RData$", ".wig", args$input)
#chrom.sizes.f <- "/home/rilla/nucleServ/wig_utils/sacCer3.chrom.sizes"
chrom.sizes.f <- paste0("/home/rilla/nucleServ/wig_utils/",
                        args$genome,
                        ".chrom.sizes")

message("writing wig file")
writeWig(splited, wigf)
message("converting to bigWig")
system(paste(tobig.bin, wigf, chrom.sizes.f, args$output))
file.remove(wigf)

###############################################################################
