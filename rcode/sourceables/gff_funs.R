#!/usr/bin/Rscript

library(IRanges)
library(GenomicRanges)

rd2df <- function (rd)
{   # Convert a RangedData into a data.frame to ease later conversion into
    # gff format.
    # It really just sets the the "space" field to "seqname", and removes the
    # field "width". Leaves everything else unchanged.
    df <- as.data.frame(rd)
    names(df) <- sub("space", "seqname", names(df))
    df[["width"]] <- NULL
    df
}

df2rd <- function (df)
{
    names(df) <- sub("seqname", "space", names(df))
    RangedData(df)
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

writeGff <- function (df, outpath)
{
    # Use this to write the output of df2gff to disk.
    write.table(df,
                file=outpath,
                quote=FALSE,
                sep="\t",
                row.names=FALSE,
                col.names=FALSE)
}

readGff <- function (fname)
{
    cols <- c("seqname",
              "source",
              "feature",
              "start",
              "end",
              "score",
              "strand",
              "frame",
              "attribute")

    df <- read.table(fname,
                     col.names=cols,
                     stringsAsFactors=FALSE)

    for (i in colnames(df)) {
        if (all(df[[i]] == ".")) {
            df[[i]] <- NULL
        }
    }
    attrs <- df$attribute
    df$attribute <- NULL

    parseRowAttrs <- function(x) {
        vals <- c()
        for (x in strsplit(x, "=")) {
            vals[[x[1]]] <- x[2]
        }
        vals
    }

    cbind(df,
          as.data.frame(do.call(rbind,
                                lapply(strsplit(attrs, ";"),
                                       parseRowAttrs))))
}
