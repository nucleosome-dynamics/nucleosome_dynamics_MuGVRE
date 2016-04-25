#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/coverage.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.RData"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################
