#!/bin/bash

bfile=$1
avdep=$2
mindep=$((avdep/3))
#mindep=13
maxdep=$((2*avdep))
nthr=$3
mkdir sfs

li=bamlist.txt
rm $li
echo $bfile >> $li
if [ ! -f $bfile.bai ]; then
samtools index $bfile
fi

out=$(basename $bfile .bam)
{ time /home/debora/src/angsd/angsd -C 50 -rf ~/space2/misti/human/strict_mask/strictmask_anccons_autosomes.regions -sites ~/space2/misti/human/strict_mask/strictmask_anccons_autosomes.sites -setMinDepth $mindep -setMaxDepth $maxdep -GL 1 -minMapQ 30 -minQ 30 -b $li -anc /space/s2/debora/misti/human/ancestral/Ancestral_states/anccons_allchr.fa -ref ~/shared/human/refs/hs37d5.fa -P $nthr -out sfs/${out} -doSaf 1 -doCounts 1 ; } 2> sfs/dosaf_${out}.time 
