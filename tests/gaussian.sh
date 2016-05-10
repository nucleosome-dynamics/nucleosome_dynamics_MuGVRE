#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_nr.gff"
reads="/home/rilla/scratch/nucler/cell_cycle_data/reads/g1.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/stiff/g1chrI.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output \
    --chr "chrI" &

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nr.gff"
reads="/home/rilla/scratch/nucler/cell_cycle_data/reads/s.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/stiff/schrI.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output \
    --chr "chrI" &

###############################################################################

wait
