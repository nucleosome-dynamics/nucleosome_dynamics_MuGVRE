#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nucleosomeDynamics.R"

input1="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep2_30m_S.bam.RData"
input2="/orozco/services/Rdata/Web/USERS/ND577a8fb9e334c/uploads/rep2_00m_G1.bam.RData"

outputGff="~/hs.gff"
plotRData="~/prd.RData"
dynRData="drd.RData"

cores="1"
equalSize="FALSE"
roundPow="5"
readSize="140"
maxDiff="70"
maxLen="140"
same_magnitude="2"
combined="TRUE"

Rscript $script \
    --input1          $input1         \
    --input2          $input2         \
    --outputGff       $outputGff      \
    --plotRData       $plotRData      \
    --dynRData        $dynRData       \
    --cores           $cores          \
    --equalSize       $equalSize      \
    --roundPow        $roundPow       \
    --readSize        $readSize       \
    --maxDiff         $maxDiff        \
    --maxLen          $maxLen         \
    --same.magnitude  $same_magnitude \
    --combined        $combined
