Making psmc:

./psmc.all.sh SAMPLE #pull 4 dependent script to cluster to make psmc


Run angsd doSaf

sbatch -A proj_1669 --wrap="./angsd.sh SAMPLE nthr" -c nthr #nthr=4 is optimal


Make 2dSFS:

sbatch -A proj_1669 --wrap"./window.2d.sfs SAMPLE1 SAMPLE2 nthr" -c
