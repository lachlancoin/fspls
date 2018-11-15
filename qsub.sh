#!/bin/bash
#PBS -S /bin/bash
#PBS -A Group-Coin
#PBS -l nodes=1:ppn=1,mem=256G,vmem=256G
#PBS -l walltime=48:00:00


cd /shares/common/groups/Group-Coin/c.zhou/fspls/fs-pls

module load R

time Rscript cmd.R > OFS.out 2> OFS.err
