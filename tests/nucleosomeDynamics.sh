#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/nucleosomeDynamics.R"
input1="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-34.RData"
input2="/orozco/services/Rdata/tmp_wd/120502_SN365_B_L002_GGM-35.RData"
plotRData="/orozco/services/Rdata/tmp_wd/nd_34_35_plot.RData"
dynRData="/orozco/services/Rdata/tmp_wd/nd_34_35.RData"
outputGff="/orozco/services/Rdata/tmp_wd/nd_34_35.gff"
#threshRData
#rep1
#rep2

Rscipt $script             \
    --input1 $input1       \
    --input2 $input2       \
    --plotRData $plotRData \
    --outputGff $outputGff \
    --dynRData $dynRData


###############################################################################
