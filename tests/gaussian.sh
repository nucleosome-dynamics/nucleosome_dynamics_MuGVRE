#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_NR.gff"
reads="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
start=1000
end=10000
chr="chrI"
output="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34_gauss.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --start $start   \
    --output $output \
    --end $end       \
    --chr $chr

###############################################################################
