#!/bin/bash

#sbatch -A proj_1669 --wrap="./make.fq.sh HGDP00666 4" -c 4

echo "MAKE FQ"


SAMPLE=$1
nthr=$2
ref_fasta=/home/av.ilina/datasets/fasta/ref_genome/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa

CRAM=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram

avdep=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' ./${SAMPLE}/${SAMPLE}.coverage.txt)
printf "average coverage is "
echo $avdep
mindep=$((avdep/3))
maxdep=$((2*avdep))


cd ${SAMPLE}
mkdir -p fq

process_chr_fq() {
    chr=$1
    local chr_name="chr${chr}"
    local cram_chr="./${chr_name}.cram"
 
    local bed="/home/av.ilina/datasets/masks/grch38/autosome_${chr_name}.bed"   

    bcftools mpileup -C50 -f "${ref_fasta}"  -T ${bed} ${cram_chr} | \
    bcftools call -c - |  \
    vcfutils.pl vcf2fq -d $mindep -D $maxdep  -Q 20| gzip > fq/${SAMPLE}.chr${chr}.fq.gz 
}



export -f process_chr_fq
export ref_fasta SAMPLE mindep maxdep nthr BED


parallel -j ${nthr} process_chr_fq ::: {1..22} 

echo "All chromosomes processed"

wait
cat fq/${SAMPLE}.chr*.fq.gz  > fq/${SAMPLE}.fq.gz && rm fq/${SAMPLE}.chr*.fq.gz && rm *.cram*
