#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nfr.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_nfr.gff"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################

script="/home/rilla/nucleServ/bin/nfr.R"
input="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35_NR.gff"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35_nfr.gff"

Rscript $script    \
    --input $input \
    --output $output

###############################################################################
