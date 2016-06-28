#!/usr/bin/env sh

###############################################################################

#script="/home/rilla/nucleServ/bin/js_plot.R"
script="/usr/people/rilla/nucleServ/bin/js_plot.R"
input="/orozco/services/Rdata/tmp_wd/g1_s_nd_plot.RData"
output="/usr/people/rilla/public_html/test_plot.html"
start="27000"
end="28000"
chr="chrI"

Rscript $script      \
    --input $input   \
    --output $output \
    --start $start   \
    --end $end       \
    --chr $chr

###############################################################################
