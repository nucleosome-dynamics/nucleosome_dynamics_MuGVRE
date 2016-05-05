#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/to_bigWig.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_cov.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_cov.bw"
genome="R64-1-1"

Rscript $script      \
    --input $input   \
    --output $output \
    --genome $genome

###############################################################################

script="/home/rilla/nucleServ/bin/to_bigWig.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_cov.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_cov.bw"
genome="R64-1-1"

Rscript $script      \
    --input $input   \
    --output $output \
    --genome $genome

###############################################################################
