#!/usr/bin/env sh

###############################################################################

script="/orozco/services/Rdata/Web/apps/nucleServ/bin/nucleosomeDynamics.R"

input1="/orozco/services/Rdata/Web/USERS/ND577a949dbd7b8/uploads/ih_methyl_V136_noMet.rep1_sorted.bam.RData"
input2="/orozco/services/Rdata/Web/USERS/ND577a949dbd7b8/uploads/ih_methyl_V138_met28h.rep1_sorted.bam.RData"
rep1="/orozco/services/Rdata/Web/USERS/ND577a949dbd7b8/uploads/ih_methyl_V137_noMet.rep2_sorted.bam.RData"
rep2="/orozco/services/Rdata/Web/USERS/ND577a949dbd7b8/uploads/ih_methyl_V139_met28h.rep2_sorted.bam.RData"

outputGff="ND_ih_methyl_V136_noMet.rep1_sorted-ih_methyl_V138_met28h.rep1_sorted.gff"
plotRData="ND_ih_methyl_V136_noMet.rep1_sorted-ih_methyl_V138_met28h.rep1_sorted_plot.RData"
dynRData="ND_ih_methyl_V136_noMet.rep1_sorted-ih_methyl_V138_met28h.rep1_sorted.RData"
threshRData="ND_ih_methyl_V136_noMet.rep1_sorted-ih_methyl_V138_met28h.rep1_sorted_thr.RData"
logf="ND_4-2.log"

cores="1"
equalSize="FALSE"
roundPow="5"
readSize="140"
maxDiff="70"
maxLen="140"
same_magnitude="2"
threshold="4"
combined="TRUE"

Rscript /orozco/services/Rdata/Web/apps/nucleServ/bin/nucleosomeDynamics.R
    --input1          $input1         \
    --input2          $input2         \
    --outputGff       $outputGff      \
    --plotRData       $plotRData      \
    --dynRData        $dynRData       \
    --threshRData     $threshRData    \
    --rep1            $rep1           \
    --rep2            $rep2           \
    --cores           $cores          \
    --equalSize       $equalSize      \
    --roundPow        $roundPow       \
    --readSize        $readSize       \
    --maxDiff         $maxDiff        \
    --maxLen          $maxLen         \
    --same.magnitude  $same_magnitude \
    --threshold       $threshold      \
    --combined        $combined       \
    >& $logf
