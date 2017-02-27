#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/coverage.R"
input="/home/rilla/scratch/jurgen_nucs/pooled.RData"
output="/home/rilla/scratch/jurgen_nucs/pooled_cov.RData"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/to_bigWig.R"
input="/home/rilla/scratch/jurgen_nucs/pooled_cov.RData"
output="/home/rilla/scratch/jurgen_nucs/pooled_cov.wig"
genome="mm9"

Rscript $script      \
    --input $input   \
    --output $output \
    --genome $genome

###############################################################################
