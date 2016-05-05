#!/usr/bin/env sh

################################################################################
#
#script="/home/rilla/nucleServ/bin/coverage.R"
#input="/orozco/services/Rdata/tmp_wd/short_34.RData"
#output="/orozco/services/Rdata/tmp_wd/short_34_cov.RData"
#
#Rscript $script    \
#    --input $input \
#    --output $output
#
################################################################################
#
#script="/home/rilla/nucleServ/bin/coverage.R"
#input="/orozco/services/Rdata/tmp_wd/short_35.RData"
#output="/orozco/services/Rdata/tmp_wd/short_35_cov.RData"
#
#Rscript $script    \
#    --input $input \
#    --output $output
#
################################################################################

script="/home/rilla/nucleServ/bin/coverage.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/reads/g1.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/g1_cov.RData"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/coverage.R"
input="/home/rilla/scratch/nucler/cell_cycle_data/reads/s.RData"
output="/home/rilla/scratch/nucler/cell_cycle_data/calcs/s_cov.RData"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################
