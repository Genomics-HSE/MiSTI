#!/bin/bash


SAMPLE=$1


#sbatch -A proj_1669 --wrap="./split.cram.sh SAMPLE 8" -c 24

nthr_per_proc=$2

CRAM=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram
CRAI=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram.crai
ref_fasta=/home/av.ilina/datasets/fasta/ref_genome/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa


avdep=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' ./${SAMPLE}/${SAMPLE}.coverage.txt)
printf "average coverage is "
echo $avdep
mindep=$((avdep/3))
maxdep=$((2*avdep))

wait


cd ${SAMPLE}


process_cram_chr() {
  chr=$1
  chr_name="chr${chr}"
  output_cram=${chr_name}.cram

  echo "Splitting $chr_name..."
  samtools view -T "${ref_fasta}" -C -h -o "$output_cram" "${CRAM}" "$chr_name" --threads "${nthr_per_proc}" && \ 
  samtools index -c "$output_cram" --threads ${nthr_per_proc}
  echo "Created $output_cram"
}

# Export variables and function to make them available to parallel
export CRAM ref_fasta nthr_per_proc
export -f process_cram_chr

# Run in parallel (adjust -j for your CPU cores)
parallel -j 3 process_cram_chr ::: {1..22} 




