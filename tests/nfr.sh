#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nfr.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_nr.gff"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_nfr.gff"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/nfr.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nr.gff"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nfr.gff"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################
