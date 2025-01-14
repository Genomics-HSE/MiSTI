#!/bin/bash
prefix=$1
genome1=$2
genome2=$3
jsfs=$4
bs_number=100
st=125

./utils/calc_time.py ./${prefix}/${genome1} ./${prefix}/${genome2}  |wc -l

echo "done"

echo ./${prefix}/${jsfs}

./utils/generateJSFS_bs.py ${bs_number} ./${prefix}/${jsfs} > ./${prefix}/bs.12.sfs


echo "done"

LLH=LLH.12

cd ${prefix}
mkdir ${LLH}
cd ../


j=3

for bs in {0..100};
do
for st in {22..30};
do
#./MiSTI.py ./${prefix}/${genome1} ./${prefix}/${genome2}  ./${prefix}/bs.sfs ${st} -bs ${bs}|grep "llh"|  awk '{print $14}' >> ./${prefix}/LLH/LLH.bs=${bs}.txt
./MiSTI.py ./${prefix}/${genome1} ./${prefix}/${genome2}  ./${prefix}/bs.12.sfs ${st} -bs ${bs}  -uf -mi 1 4 ${st} ${j} 1 --cpfit|grep "llh"|tr -d "][" |  awk '{print $17,$14,$6} ' >> ./${prefix}/${LLH}/LLH.bs=${bs}.txt

#./MiSTI.py ./${prefix}/${genome1} ./${prefix}/${genome2}  ./${prefix}/bs.sfs ${st} -bs ${bs}  -uf -mi 1 4 ${st} ${j} 1 -mi 2 4 ${st} ${j} 1  --cpfit|grep "llh"|tr -d "][" |  awk '{print $18,$14, $15,$6} ' >> ./${prefix}/${LLH}/LLH.bs=${bs}.txt
done
done
