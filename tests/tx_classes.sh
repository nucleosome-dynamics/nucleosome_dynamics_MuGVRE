#!/usr/bin/env sh

################################################################################
#
#script="/home/rilla/nucleServ/bin/tx_classes.R"
#calls="/orozco/services/Rdata/tmp_wd/short_35_NR.gff" #coverage="/orozco/services/Rdata/tmp_wd/short_35_cov.RData"
#output="/orozco/services/Rdata/tmp_wd/short_test_tx.gff"
#genome="R64-1-1"
#
#Rscript $script          \
#    --calls $calls       \
#    --coverage $coverage \
#    --genome $genome     \
#    --output $output     \
#    --cores 20
#
################################################################################

#script="/home/rilla/nucleServ/bin/tx_classes.R"
#calls="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/calcs/s_nr.gff"
#coverage="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_cov.RData"
#output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_tss_new.gff"
#genome="R64-1-1"

script="/home/rilla/nucleServ/bin/tx_classes.R"
calls="/orozco/scratch/xesh0/orozco/rilla/nucler/cell_cycle_data/calcs/g1_nr.gff"
coverage="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_cov.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_tss_new.gff"
genome="R64-1-1"

Rscript $script          \
    --calls $calls       \
    --coverage $coverage \
    --genome $genome     \
    --output $output     \
    --cores 20

###############################################################################
