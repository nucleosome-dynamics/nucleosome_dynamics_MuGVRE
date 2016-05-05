#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_nr.gff"
reads="/home/rilla/scratch/nucler/cell_cycle_data/reads/g1.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/reads/gau_g1.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nr.gff"
reads="/home/rilla/scratch/nucler/cell_cycle_data/reads/s.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/reads/gau_s.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output

###############################################################################
