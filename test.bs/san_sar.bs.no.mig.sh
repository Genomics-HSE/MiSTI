#!/bin/bash
prefix=$1
genome1=$2
genome2=$3
jsfs=$4
bs_number=1000
st=125

./utils/calc_time.py ./${prefix}/${genome1} ./${prefix}/${genome2}  |wc -l

echo "done"

echo ./${prefix}/${jsfs}

./utils/generateJSFS_bs.py ${bs_number} ./${prefix}/${jsfs} > ./${prefix}/2.bs.no.mig.sfs


echo "done"

LLH=LLH.no.mig.2

cd ${prefix}
mkdir ${LLH}
cd ../


j=3

for bs in {0..1000};
do
for st in {15..23};
do
./MiSTI.py ./${prefix}/${genome1} ./${prefix}/${genome2}  ./${prefix}/2.bs.no.mig.sfs ${st} -bs ${bs}|grep "llh"|  awk '{print $14,$6}' >> ./${prefix}/${LLH}/LLH.bs=${bs}.txt

done
done
