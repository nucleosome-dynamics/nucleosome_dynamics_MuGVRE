#!/usr/bin/env sh

################################################################################
#
#script="/home/rilla/nucleServ/bin/readBAM.R"
#type="paired"
#input="/orozco/services/Rdata/in_bams/120502_SN365_B_L002_GGM-34.bam"
#output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
#
#Rscript $script    \
#    --type $type   \
#    --input $input \
#    --output $output
#
################################################################################

script="/home/rilla/nucleServ/bin/readBAM.R"
type="paired"
input="/orozco/services/Rdata/in_bams/120502_SN365_B_L002_GGM-35.bam"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35.RData"

Rscript $script    \
    --type $type   \
    --input $input \
    --output $output

###############################################################################
