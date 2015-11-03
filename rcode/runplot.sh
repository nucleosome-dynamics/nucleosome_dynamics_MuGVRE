#!/bin/sh

dir="/orozco/services/R-data/tmp_wd"
name="120502_SN365_B_L002_GGM-34_120502_SN365_B_L002_GGM-35"
start=1000
end=2000
chr="chrI"

Rscript plot.R \
    --start $start \
    --end $end \
    --chr $chr \
    --input ${dir}/${name}.RData \
    --output ${dir}/${name}_${start}-${end}_${chr}.png
