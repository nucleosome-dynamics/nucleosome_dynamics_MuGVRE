#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/readBAM.R"
type="paired"
input="/orozco/services/Rdata/in_bams/short_34.bam"
output="/orozco/services/Rdata/tmp_wd/short_34.RData"

Rscript $script    \
    --type $type   \
    --input $input \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/readBAM.R"
type="paired"
input="/orozco/services/Rdata/in_bams/short_35.bam"
output="/orozco/services/Rdata/tmp_wd/short_35.RData"

Rscript $script    \
    --type $type   \
    --input $input \
    --output $output

###############################################################################
