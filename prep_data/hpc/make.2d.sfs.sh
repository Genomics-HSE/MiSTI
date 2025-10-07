#!/bin/bash

#./make.2d.sfs.sh nthr sample1 sample2

nthr=$1
sample1=$2
sample2=$3
pop1=$4
pop2=$5


mkdir -p ./2dsfs
realSFS  ./${sample1}/${sample1}.sfs/${sample1}.all.saf.idx   ./${sample2}/${sample2}.sfs/${sample2}.all.saf.idx -cores ${nthr}  -nSites 2500000  > ./2dsfs/${sample1}_${sample2}.new.2dsfs             

wait

./utils/ANGSDSFS.py ./2dsfs/${sample1}_${sample2}.new.2dsfs ${pop1} ${pop2} > ./2dsfs/${sample1}_${sample2}.new.sfs 
