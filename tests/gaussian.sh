#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/gauss_fit.R"
calls="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/run031/NR_rep1_00m_G1.gff"
reads="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep1_00m_G1.bam.RData"
output="/home/rilla/stf.gff"

Rscript $script      \
    --calls $calls   \
    --reads $reads   \
    --output $output

###############################################################################
