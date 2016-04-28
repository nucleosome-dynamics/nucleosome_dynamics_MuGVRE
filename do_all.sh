#!/usr/bin/env sh

wd="/home/rilla/scratch/nucler/cell_cycle_data/workbench"
bin_dir="/home/rilla/scratch/nucler/nucleServ/bin"

genome="R64-1-1"

rs1="${wd}/g1_reads.RData"
cov1="${wd}/g1_cov.RData"
nr1="${wd}/g1_calls.gff"
tx_class1="${wd}/g1_txclasses.gff"
nfr1="${wd}/g1_nfr.gff"
per_bw1="${wd}/g1_per.bw"
per_gff1="${wd}/g1_per.gff"
gauss1="${wd}/g1_gauss.gff"

rs2="${wd}/s_reads.RData"
cov2="${wd}/s_cov.RData"
tx_class2="${wd}/s_txclasses.gff"
nr2="${wd}/s_calls.gff"
nfr2="${wd}/s_nfr.gff"
per_bw2="${wd}/s_per.bw"
per_gff2="${wd}/s_per.gff"
gauss2="${wd}/s_gauss.gff"

nd="${wd}/g1_s_nd.gff"


echo "========================================================================"
echo "coverage"
echo "========================================================================"

script="${bin_dir}/coverage.R"

Rscript $script  \
    --input $rs1 \
    --output $cov1

Rscript $script  \
    --input $rs2 \
    --output $cov2


echo "========================================================================"
echo "nucleR"
echo "========================================================================"

script="${bin_dir}/nucleR.R"

Rscript $script  \
    --input $rs1 \
    --output $nr1

Rscript $script  \
    --input $rs2 \
    --output $nr2


echo "========================================================================"
echo "TX Classes"
echo "========================================================================"

script="${bin_dir}/tx_classes.R"

Rscript $script      \
    --calls $nr1     \
    --coverage $cov1 \
    --genome $genome \
    --output $tx_class1

Rscript $script      \
    --calls $nr2     \
    --coverage $cov2 \
    --genome $genome \
    --output $tx_class2


echo "========================================================================"
echo "NFR"
echo "========================================================================"

script="${bin_dir}/nfr.R"

Rscript $script  \
    --input $nr1 \
    --output $nfr1

Rscript $script  \
    --input $nr2 \
    --output $nfr2


echo "========================================================================"
echo "periodicity"
echo "========================================================================"

script="${bin_dir}/peiodicity.R"

Rscript $script          \
    --calls $nr1         \
    --genome $genome     \
    --bwOutput $per_bw1  \
    --gffOutput $gff_bw1 \
    --coverage $cov1

Rscript $script          \
    --calls $nr2         \
    --genome $genome     \
    --bwOutput $per_bw2  \
    --gffOutput $gff_bw2 \
    --coverage $cov2


echo "========================================================================"
echo "Gaussian fitting"
echo "========================================================================"

script="${bin_dir}/gauss_fit.R"

Rscript $script  \
    --calls $nr1 \
    --reads $rs1 \
    --output $gauss1

Rscript $script  \
    --calls $nr2 \
    --reads $rs2 \
    --output $gauss2


echo "========================================================================"
echo "Nucleosome Dynamics"
echo "========================================================================"

script="${bin_dir}/nucleosomeDynamics.R"

Rscript $script   \
    --input1 $rs1 \
    --input2 $rs2 \
    --outputGff $nd
