#!/bin/bash


#./coverage.sh sample nthr
SAMPLE=$1 # sample name
CRAM=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram
CRAI=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram.crai


nthr=$2
mkdir -p ${SAMPLE}

ref_fasta=/home/av.ilina/datasets/fasta/ref_genome/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa

if [ ! -f ${CRAI} ]; then
samtools index ${CRAM}
fi

# calculate coverage from the cram file
samtools view -h ${CRAM} -T ${ref_fasta} --threads ${nthr}| grep -m 22 "@SQ" |sed 's/:/\t/g'|awk '{print $3"\t"0"\t"$5}' > ./${SAMPLE}/${SAMPLE}.autosomes.bed
samtools bedcov ./${SAMPLE}/${SAMPLE}.autosomes.bed  $CRAM |grep -vP "\t0$" >> ./${SAMPLE}/${SAMPLE}.coverage.txt

AVCOV=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' ./${SAMPLE}/${SAMPLE}.coverage.txt)
printf "average coverage is"
echo $AVCOV
