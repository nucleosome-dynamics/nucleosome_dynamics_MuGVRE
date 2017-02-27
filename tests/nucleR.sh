#!/usr/bin/env sh

################################################################################
#
#script="/home/rilla/nucleServ/bin/nucleR.R"
#input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
#output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
#
#Rscript $script    \
#    --input $input \
#    --output $output
#
################################################################################

###############################################################################

script="/home/rilla/nucleServ/bin/nucleR.R"
#input="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/reads/g1.RData"
#output="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/calcs/g1_nr.gff"
input="/home/rilla/scratch/jurgen_nucs/pooled.RData"
output="/home/rilla/scratch/jurgen_nucs/pooled_calls.gff"

Rscript $script    \
    --input $input \
    --output $output \
    --cores 20

###############################################################################
#
#script="/home/rilla/nucleServ/bin/nucleR.R"
#input="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/reads/s.RData"
#output="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/calcs/s_nr.gff"
#
#Rscript $script    \
#    --input $input \
#    --output $output \
#    --cores 20
#
################################################################################
