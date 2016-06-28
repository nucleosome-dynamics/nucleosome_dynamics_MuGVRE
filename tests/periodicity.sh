#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/periodicity.R"
calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_nr.gff"
coverage="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_cov.RData"
genome="R64-1-1"
bwOutput="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_period.bw"
gffOutput="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_period.gff"

Rscript $script            \
    --calls $calls         \
    --coverage $coverage   \
    --genome $genome       \
    --bwOutput $bwOutput   \
    --gffOutput $gffOutput

################################################################################
#
#script="/home/rilla/nucleServ/bin/periodicity.R"
#calls="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_nr.gff"
#coverage="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_cov.RData"
#genome="R64-1-1"
#bwOutput="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_period.bw"
#gffOutput="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_period.gff"
#
#Rscript $script            \
#    --calls $calls         \
#    --coverage $coverage   \
#    --genome $genome       \
#    --bwOutput $bwOutput   \
#    --gffOutput $gffOutput
#
################################################################################
