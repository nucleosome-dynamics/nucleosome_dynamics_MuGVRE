#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/from_bigWig.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.bw"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov_new.RData"

Rscript $script      \
    --input $input   \
    --output $output

###############################################################################
