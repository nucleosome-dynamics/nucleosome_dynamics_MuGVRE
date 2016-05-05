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



###############################################################################

script="/home/rilla/nucleServ/bin/coverage.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35.RData"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35_cov.RData"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################
