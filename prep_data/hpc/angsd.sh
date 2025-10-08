
#!/bin/bash


 
#calculate join sfs
#./angsdsfs.sh AVDEP nTHR SAMPLE

avdep=$1
nthr=$2
SAMPLE=$3


if [[ "$avdep" =~ ^[0-9]+$ ]]; then
    echo "AVDEPTH: $avdep"
else
    avdep=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' ${SAMPLE}.coverage.txt)
        printf "average coverage for ${SAMPLE1} is "
        echo $avdep
fi

mindep=$((avdep/3))
maxdep=$((2*avdep))

anc_fasta=/home/av.ilina/datasets/Ancestral/homo_sapiens_ancestor_GRCh38/ancestral.grch38.autosomes.fa
ref_fasta=/home/av.ilina/datasets/fasta/ref_genome/grch38/GRCh38_full_analysis_set_plus_decoy_hla.fa
CRAM=/home/av.ilina/datasets/BAMs/${SAMPLE}.cram
strict_mask_sites=/home/av.ilina/datasets/masks/grch38/autosome.sites.strict.mask.bed


mkdir -p ${SAMPLE}/${SAMPLE}.sfs


angsd -P $nthr -i ${CRAM} -out ./${SAMPLE}/${SAMPLE}.sfs/${SAMPLE}.all -sites ${strict_mask_sites} \
 -C 50 \
 -minMapQ 20 -minQ 30 \
 -setMinDepth $mindep -setMaxDepth $maxdep \
 -GL 1 -ref ${ref_fasta}  -anc ${anc_fasta} \
 -doSaf 1 \
 -doCounts 1





