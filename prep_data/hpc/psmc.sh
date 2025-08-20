#!/bin/bash


#sbatch -A proj_1669 --wrap="./psmc.sh HGDP00666" 
SAMPLE=$1

cd ${SAMPLE}
mkdir -p psmcIN
mkdir -p psmcOUT

echo "Running fq2psmcfa"
{ time fq2psmcfa -q20 fq/${SAMPLE}.fq.gz > psmcIN/${SAMPLE}.psmcfa ; } 2> psmcIN/fq2psmc_${SAMPLE}.time
wait

echo " we finishes psmcfa..."
{ time psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o psmcOUT/${SAMPLE}.psmc psmcIN/${SAMPLE}.psmcfa ; } 2> psmcOUT/psmc_${SAMPLE}.time
wait




# PSMC plot
psmc_plot.pl -n25 -u 1.25e-8 -g1 -x1000 -X1000000 -L psmcOUT/plot_${SAMPLE} psmcOUT/${SAMPLE}.psmc



