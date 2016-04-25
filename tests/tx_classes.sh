#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/tx_classes.R"
calls="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
coverage="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.RData"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_txclasses.gff"
genome="R64-1-1"

Rscript $script          \
    --calls $calls       \
    --coverage $coverage \
    --genome $genome     \
    --output $output
