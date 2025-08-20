#!/bin/bash

#sbatch -A proj_1669 --wrap="./make.fq.sh HGDP00666 4" -c 4
SAMPLE=$1
nthr=$2
ref_fasta=/home/av.ilina/datasets/fasta/ref_genome/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa


avdep=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' ./${SAMPLE}/${SAMPLE}.coverage.txt)
printf "average coverage is "
echo $avdep
mindep=$((avdep/3))
maxdep=$((2*avdep))
wait


export mindep maxdep ref_fasta
cd ${SAMPLE}

mkdir -p fq

process_chr() {
    local chr=$1
    local BED="/home/av.ilina/datasets/masks/grch38/autosome_chr${chr}.bed"
    local CRAM="/home/av.ilina/MiSTI/misti_clone/${SAMPLE}/chr21.cram"
    
    
    # Variant calling pipeline
    bcftools mpileup -C50 -f "${ref_fasta}" -R "${BED}" "${CRAM}" --threads 2 -Ou | \
    bcftools call -mv --threads 2 |bcftools filter -i 'QUAL>=20'  -Oz -o "${SAMPLE}.chr${chr}.vcf.gz"
    
    # Index and generate consensus
    bcftools index "${SAMPLE}.chr${chr}.vcf.gz" --threads 2  && \
    bcftools consensus -m "${BED}" -f "${ref_fasta}" -i "DP>${mindep} && DP<${maxdep}" \
        "${SAMPLE}.chr${chr}.vcf.gz" | gzip > "fq/${SAMPLE}_chr${chr}.fq.gz"
    
    echo "Finished chromosome ${chr} at $(date)"
}

# Export variables and function
export mindep maxdep SAMPLE
export -f process_chr

# Process chromosomes 21 and 22 in parallel
echo "Starting parallel processing at $(date)"
parallel -j 2 --linebuffer --tag "process_chr {}" ::: 21 22
echo "All chromosomes processed at $(date)" 

wait
echo "All chromosomes processed"




wait
cat fq/${SAMPLE}_chr*.fq.gz  > fq/${SAMPLE}.fq.gz
#rm fq/${SAMPLE}_chr*.fq.gz

