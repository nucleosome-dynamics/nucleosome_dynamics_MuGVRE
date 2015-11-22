#!/bin/bash

interp="/usr/bin/Rscript"
script="/home/rilla/nucleServ/rcode/coverage.R"
wd="/orozco/services/Rdata/tmp_wd"
inf="120502_SN365_B_L002_GGM-34.RData"
outf="120502_SN365_B_L002_GGM-34_cov.RData"
cores=1
type="paired"

$interp $script            \
    --input ${wd}/${inf}   \
    --output ${wd}/${outf} \
    --cores $cores         \
    --type $type
