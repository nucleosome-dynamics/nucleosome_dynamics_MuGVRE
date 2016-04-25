#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/periodicity.R"
calls="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
coverage="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_cov.RData"
genome="R64-1-1"
bwOutput="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_period.bw"
gffOutput="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_period.gff"

Rscript $script            \
    --calls $calls         \
    --coverage $coverage   \
    --genome $genome       \
    --bwOutput $bwOutput   \
    --gffOutput $gffOutput

###############################################################################
