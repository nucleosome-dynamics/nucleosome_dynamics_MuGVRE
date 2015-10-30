#!/bin/sh

/usr/bin/Rscript /home/rilla/nucleServ/rcode/nucleosomeDynamics.R \
    --input1 /orozco/services/R-data/tmp_wd/120502_SN365_B_L002_GGM-34.RData \
    --input2 /orozco/services/R-data/tmp_wd/120502_SN365_B_L002_GGM-35.RData \
    --outputGff /orozco/services/R-data/tmp_wd/120502_SN365_B_L002_GGM-34_120502_SN365_B_L002_GGM-35.gff \
    --outputRData /orozco/services/R-data/tmp_wd/120502_SN365_B_L002_GGM-34_120502_SN365_B_L002_GGM-35.RData \
    --cores 1 \
    --same.magnitude 2 \
    --equalSize FALSE \
    --roundPow 7 \
    --combined F \
    --maxLen 170 \
    --maxDiff 70 \
    --readSize 140 \
    --threshold 60%
