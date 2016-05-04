#!/usr/bin/env sh

###############################################################################

script="/home/rilla/nucleServ/bin/js_plot.R"
input="/orozco/services/Rdata/tmp_wd/nd_34_35_plot.RData"
output="/orozco/services/Rdata/tmp_wd/nd_34_35_plot.html"
start="5000"
end="6000"
chr="chrI"

Rscript $script      \
    --input $input   \
    --output $output \
    --start $start   \
    --end $end       \
    --chr $chr

###############################################################################
