#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/to_bigWig.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.RData"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.bw"
genome="R64-1-1"

Rscript $script      \
    --input $input   \
    --output $output \
    --genome $genome

###############################################################################
