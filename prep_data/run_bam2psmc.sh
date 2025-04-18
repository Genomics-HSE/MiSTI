#!/bin/bash

bfile=$1
avdep=$2
mindep=$((avdep/3))
maxdep=$((2*avdep))
ncores=$3

mkdir fq
mkdir psmcIN
mkdir psmcOUT

out=$(basename $bfile .bam)
for chr in {1..22}
do
  echo "{ time samtools mpileup -C50 -u -v -f /space/s2/sharedResources/human/refs/hs37d5.fa --positions /space/s2/debora/misti/human/strict_mask/strictmask_anccons_chr${chr}.bed $bfile 2> fq/${out}_stderr.txt  | bcftools call -c - | vcfutils.pl vcf2fq -d $mindep -D $maxdep | gzip > fq/${out}_chr${chr}.fq.gz ; } 2> fq/bam2fq_${out}_${chr}.time "
done | parallel -j $ncores
echo "Waiting for bam2fq to finish"
wait
cat fq/${out}_chr*.fq.gz  > fq/${out}.fq.gz
rm fq/${out}_chr*.fq.gz

echo "Running fq2psmcfa"
{ time fq2psmcfa -q30 fq/$out.fq.gz > psmcIN/$out.psmcfa ; } 2> psmcIN/fq2psmc_${out}.time
wait
{ time psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o psmcOUT/$out.psmc psmcIN/$out.psmcfa ; } 2> psmcOUT/psmc_${out}.time
wait
