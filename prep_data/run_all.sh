# script to run PSMC and joint SFS from human bam files

BAM=$1 # path to the bam file
PREF=$(basename $BAM .bam)

if [ ! -f ${BAM}.bai ]; then
samtools index $BAM
fi
# calculate coverage from the bam file
samtools view -h $BAM | grep -m 22 "@SQ" |sed 's/:/\t/g'|awk '{print $3"\t"0"\t"$5}' > autosomes.bed
samtools bedcov autosomes.bed $BAM |grep -vP "\t0$" >> coverage.txt

AVCOV=$(awk '{total+=$3; spl1+=$4} END {print int(spl1/total)}' coverage.txt)
printf "average coverage is"
echo $AVCOV

# PSMC with strict filter
nice run_bam2psmc.sh $BAM $AVCOV 8
# PSMC plot
psmc_plot.pl -n25 -u 1.25e-8 -g1 -x1000 -X1000000 -L psmcOUT/plot_${PREF} psmcOUT/${PREF}.psmc

# SFS with strict filter
nice run_angsdsfs.sh $BAM $AVCOV 2

