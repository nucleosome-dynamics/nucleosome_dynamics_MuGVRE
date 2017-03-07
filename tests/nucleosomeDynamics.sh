#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nucleosomeDynamics.R"

#input1="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep2_30m_S.bam.RData"
#input2="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep2_00m_G1.bam.RData"

input1="/home/rilla/scratch/nucler/cell_cycle_data/aligned/RData/rep2_00m_G1.RData"
input2="/home/rilla/scratch/nucler/cell_cycle_data/aligned/RData/rep2_30m_S.RData"

genome="R64-1-1"

outputGff="~/ND_new.gff"
outputBigWig="~/pvals.bw"
plotRData="~/prd.RData"

equalSize="FALSE"
readSize="140"
maxDiff="70"
maxLen="140"

cores="1"

Rscript $script \
    --input1          $input1         \
    --input2          $input2         \
    --genome          $genome         \
    --outputGff       $outputGff      \
    --outputBigWig    $outputBigWig   \
    --plotRData       $plotRData      \
    --equalSize       $equalSize      \
    --readSize        $readSize       \
    --maxDiff         $maxDiff        \
    --maxLen          $maxLen         \
    --cores           $cores
