#!/bin/bash

#./make.2d.sfs.sh nthr sample1 sample2

nthr=$1
sample1=$2
sample2=$3

realSFS  ./${sample1}.sfs/${sample1}.all.saf.idx   ./${sample2}.sfs/${sample2}.all.saf.idx -cores ${nthr}  -win 2500000 -step 2500000  > ${sample1}_${sample2}.2dsfs             
 
