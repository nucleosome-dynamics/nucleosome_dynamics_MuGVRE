#!/usr/bin/Rscript

## Imports ####################################################################

library(getopt)
library(IRanges)

## Parameters and Arguments ###################################################

towig.bin <- "/home/rilla/nucleServ/wig_utils/bigWigToWig"

defaults <- list(type           = NULL,
                 fdrOverAmp     = 0.05,
                 components     = 1,
                 fragmentLen    = NULL,
                 trim           = 50,
                 pcKeepComp     = 0.02,
                 width          = 125,
                 threshold      = "35%",
                 dyad.length    = NULL,
                 min.overlap    = NULL,
                 score_w.thresh = 0.6,
                 score_h.thresh = 0.4)

spec <- matrix(c("input",       "i", 1, "character",
                 "output",      "o", 1, "character"),
               byrow=TRUE,
               ncol=4)
args <- getopt(spec)

## Convert big wig to wig #####################################################

wigf <- sub(".bw$", ".wig", args$input)

message("running bigWigToWig")
system(paste(towig.bin, args$input, wigf))
message("loading wig file")
lines <- readLines(wigf)
message("removing temporary wig file from disk")
file.remove(wigf)

## Parse the wig file into a list of Rle objects ##############################

message("parsing lines from wig file")
sep.idxs <- grep("^fixedStep", lines)

vals <-  mapply(function(i, j) as.numeric(lines[i:j]),
                sep.idxs + 1,
                c(sep.idxs[-1] - 1, length(lines)),
                SIMPLIFY=FALSE)

ids <- lines[sep.idxs]

xs <- unlist(strsplit(ids, " "))

getVal <- function(xs, a)
    sub(paste0(a, "="), "", grep(a, xs, value=TRUE))

chrs <- getVal(xs, "chrom")
pos <- as.numeric(getVal(xs, "start"))

chroms <- unique(chrs)
cover <- lapply(
    chroms,
    function(chr) {
        chr.loc <- chrs == chr
        chr.pos <- pos[chr.loc]
        chr.vals <- vals[chr.loc]
        max.pos <- max(chr.pos)
        tail.size <- length(chr.vals[[which(chr.pos == max.pos)]]) - 1
        s <- max.pos + tail.size
        x <- rep(0, s)
        for (i in seq_along(chr.vals)) {
            from <- chr.pos[i]
            to <- from + length(chr.vals[[i]]) - 1
            x[from:to] <- chr.vals[[i]]
        }
        Rle(x)
    }
)
names(cover) <- chroms

## Save as an RData ###########################################################

message("stroing output as RData")
save(cover, file=args[["output"]])

###############################################################################
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
