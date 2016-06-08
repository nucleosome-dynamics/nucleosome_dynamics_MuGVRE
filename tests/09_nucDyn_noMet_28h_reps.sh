#!/bin/bash
######################################## PATHS
#cd /orozco/scratch/xesh1/orozco/dbuitrago/2014_isabelle_methylation/nucleR2

script="/home/rilla/nucleServ/bin/nucleosomeDynamics.R"

READS_FOLDER="/orozco/scratch/xesh1/orozco/dbuitrago/2014_isabelle_methylation/nucleR2/02_reads"
OUTPUT_FOLDER="/home/rilla/scratch/nucler/threshs"

mkdir -p $OUTPUT_FOLDER



###############################################################################
# Non-meth vs meth

input1=$READS_FOLDER/ih_methyl_V136_noMet.rep1.RData
input2=$READS_FOLDER/ih_methyl_V138_met28h.rep1.RData
plotRData=$OUTPUT_FOLDER/noMet1_met28h1_nd_plot.RData
outputGff=$OUTPUT_FOLDER/noMet1_met28h1_nd.gff
dynRData=$OUTPUT_FOLDER/noMet1_met28h1_nd.RData
rep1=$READS_FOLDER/ih_methyl_V137_noMet.rep2.RData
rep2=$READS_FOLDER/ih_methyl_V139_met28h.rep2.RData

Rscript $script \
--input1 $input1 \
--input2 $input2 \
--plotRData $plotRData \
--outputGff $outputGff \
--dynRData $dynRData \
--rep1 $rep1 \
--rep2 $rep2 \
--cores 10 \
--threshRData ${OUTPUT_FOLDER}/thresh.RData

###############################################################################



