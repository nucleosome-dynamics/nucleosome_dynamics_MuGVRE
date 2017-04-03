#!/usr/bin/env sh

/opt/R-3.1.2/bin/Rscript /home/rilla/nucleServ/bin/periodicity.R \
    --genes /orozco/services/Rdata/MuG/MuG_public/refGenomes/Saccharomyces_cerevisiae/R64-1-1/genes.gff \
    --reads /home/rilla/scratch/test_data/M_chrII.RData \
    --gffOutput /home/rilla/scratch/test_data/my_test/P_M_chrII.gff \
    --calls /home/rilla/scratch/test_data/my_test/NR_M_chrII.gff \
    --bwOutput /home/rilla/scratch/test_data/my_test/P_M_chrII.bw \
    --type paired \
    --chrom_sizes /orozco/services/Rdata/MuG/MuG_public/refGenomes/Saccharomyces_cerevisiae/R64-1-1/R64-1-1.fa.chrom.sizes
