#!/usr/bin/Rscript

library(IRanges)
library(GenomicRanges)

# Miscelanious helper functions

rd2df <- function(rd)
{   # Convert a RangedData into a data.frame to ease later conversion into
    # gff format.
    # It really just sets the the "space" field to "seqname", and removes the
    # field "width". Leaves everything else unchanged.
    df <- as.data.frame(rd)
    names(df) <- sub("space", "seqname", names(df))
    df[["width"]] <- NULL
    df
}

df2gff <- function (df, ...)
{   # Convert a data.frame into a form that is easily saved in disk as a gff.
    # It tries to match gff fields to the row names of the data.frame.
    # gff fields not found in the data.frame will be set to `.` and fields in
    # the data.frame not found in the gff fields will be stored accordingly in
    # the `attribute` field, separated by `;`.
    # It also accepts optional keyword arguments to specify gff fields that
    # are not present in the data.frame.
    kwargs <- list(...)
    fields <- c("seqname", "source", "feature", "start", "end", "score",
                "strand", "frame")
    out.df <- data.frame(matrix(nrow=nrow(df),
                                ncol=length(fields) + 1,
                                dimnames=list(c(),
                                              c(fields, "attribute"))))
    for (f in fields) {
        if (f %in% colnames(df)) {
            out.df[[f]] <- as.vector(df[[f]])
        } else if (f %in% names(kwargs)) {
            out.df[[f]] <- kwargs[[f]]
        } else {
            out.df[[f]] <- "."
        }
    }
    nonfield.columns <- colnames(df)[!colnames(df) %in% fields]
    attrVal <- function(i, df) sprintf("%s=%s", i, df[[i]])
    out.df[["attribute"]] <- do.call(paste,
                                     c(lapply(nonfield.columns,
                                              attrVal,
                                              df),
                                       sep=";"))
    out.df
}

writeGff <- function(df, outpath) {
    # Use this to write the output of df2gff to disk.
    write.table(df,
                file=outpath,
                quote=FALSE,
                sep="\t",
                row.names=FALSE,
                col.names=FALSE)
}

vectorizedAll <- function(...)
    # Helper function that behaves as a vectorized version of the function
    # `all`
    Reduce(`&`, list(...))

getParams <- function(defaults, inf)
{   # Given a list of default parameters and an input file for optional
    # parameters, return a list of parameters that will be used.

    parseConfig <- function(f)
    {   # Parse a tsv file containing the parameters into a named list.
        # The first column should be the parameter names and the second their
        # value.
        # Some values should be parsed as strings and others as numbers, do
        # that accordingly.
        df <- read.table(f, stringsAsFactors=FALSE)
        ls <- lapply(
            as.list(df$V2),
            function(x) {
                n <- suppressWarnings(as.numeric(x))
                if (is.na(n)) {
                    x
                } else {
                    n
                }
            }
        )
        names(ls) <- df$V1
        ls
    }
    in.params <- parseConfig(inf)
    for (p in names(in.params)) {
        defaults[[p]] <- in.params[[p]]
    }
    defaults
}

getType <- function(score_w, score_h, thresh_w, thresh_h)
    # Given a vector of score w, a vector of score h, a threshold for w and a
    # threshold for h, return a vector that tells whether those scores are
    # classified as `well-positioned`, or `fuzzy`.
    ifelse(`&`(score_w > thresh_w,
               score_h > thresh_h),
           "W", "F")

checkInF <- function(f)
{   # Some checks on the input file specified
    if (is.null(f)) {
        stop("An input file must be specified")
    }
    if (!grepl("\\.bam$", f)) {
        stop("Input file must be in BAM format")
    }
    if (!file.exists(f)) {
        stop("Specified input file doesn't exist")
    }
}

checkType <- function(t)
{   # Some checks on the sequencing type specified
    if (is.null(t)) {
        stop("A type has to be specified")
    }
    if (!t %in% c("single", "paired")) {
        stop("Type must be either `single` or `paired`")
    }
}

checkOutF <- function(f)
{   # Some checks on the output file specified
    if (is.null(f)) {
        stop("An output file must be specified")
    }
}
