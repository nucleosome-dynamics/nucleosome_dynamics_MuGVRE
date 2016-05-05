#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nucleosomeDynamics.R"
input1="/home/rilla/scratch/nucler/cell_cycle_data/reads/g1.RData"
input2="/home/rilla/scratch/nucler/cell_cycle_data/reads/s.RData"
plotRData="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_s_nd_plot.RData"
outputGff="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_s_nd.gff"

Rscript $script            \
    --input1 $input1       \
    --input2 $input2       \
    --plotRData $plotRData \
    --outputGff $outputGff \
    --threshold "80%"


###############################################################################
