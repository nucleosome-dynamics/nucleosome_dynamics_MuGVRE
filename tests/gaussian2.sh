#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nr.gff"
reads="/home/rilla/scratch/nucler/cell_cycle_data/reads/s.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/stiff/schrVIII.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output \
    --chr "chrVIII"

###############################################################################
