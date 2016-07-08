#!/usr/bin/env sh

###############################################################################

#script="/home/rilla/nucleServ/bin/js_plot.R"
script="/orozco/services/Rdata/Web/apps/nucleServ/bin/js_plot.R"
#input="/orozco/services/Rdata/tmp_wd/g1_s_nd_plot.RData"
input="/orozco/services/Rdata/Web/USERS/ND5721b9071a593/cellcycle/ND_cellcycle_G1-cellcycle_S_plot.RData"
output="/usr/people/rilla/public_html/test_plot2.html"
start="27500"
end="28500"
chr="chrI"

Rscript $script      \
    --input $input   \
    --output $output \
    --start $start   \
    --end $end       \
    --chr $chr

###############################################################################
