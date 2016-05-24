#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/readBAM.R"
type="paired"
input="/orozco/projects/yeast_methylation/round3/nucleR2/01_bam/ih_methyl_V136_noMet.rep1.bam"
output="/home/rilla/scratch/nucler/tests/ih_methyl_V136_noMet.rep1.RData"

Rscript $script    \
    --type $type   \
    --input $input \
    --output $output

###############################################################################
