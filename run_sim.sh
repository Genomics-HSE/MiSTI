#!/bin/bash

#Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 DIRECTORY \"ms arguments\"" >&2
  exit 1
fi

PSMC_PATH=/home/vladimir/psmc_project/psmc
MSHOT_PATH=/home/vladimir/psmc_project/foreign/msHOT-lite

DIR=$1
MSARG=$2

CLEAN=0
if [ "$#" -gt 2 ]; then
	if [ "$3" = "-c" ]; then
		CLEAN=1
	fi
fi

if [ -d "$DIR" ]; then
  echo "Cannot create the directory. Directory exists."
  exit 1
fi
`mkdir $DIR`
if [ ! -d "$DIR" ]; then
  echo "Failed to create the directory $DIR"
  exit 1
fi
$MSHOT_PATH/msHOT-lite $MSARG | gzip > $DIR/sim.ms.gz
./utils/MSSPLIT.py <(gunzip -c $DIR/sim.ms.gz) $DIR
gzip $DIR/ms2g1.ms
gzip $DIR/ms2g2.ms
$PSMC_PATH/utils/ms2psmcfa.pl <(gunzip -c $DIR/ms2g1.ms.gz) | gzip > $DIR/ms2g1.psmc.fa.gz
$PSMC_PATH/utils/ms2psmcfa.pl <(gunzip -c $DIR/ms2g2.ms.gz) | gzip > $DIR/ms2g2.psmc.fa.gz
#`/home/vladimir/psmc_project/psmc/psmc -p 1*4+25*2+1*4+1*6 $DIR/ms2g1.psmc.fa > $DIR/ms2g1.psm`
#`/home/vladimir/psmc_project/psmc/psmc -p 1*4+25*2+1*4+1*6 $DIR/ms2g2.psmc.fa > $DIR/ms2g2.psm`
parallel $PSMC_PATH/psmc "-p 1*4+25*2+1*4+1*6 <(gunzip -c $DIR/ms2g{}.psmc.fa.gz) > $DIR/ms2g{}.psmc" ::: 1 2
#parallel echo "-p 1*4+25*2+1*4+1*6 $DIR/ms2g{}.psmc.fa" ::: 1 2 > $DIR/ms2g{}.psmc
$PSMC_PATH/utils/psmc_plot.pl -n30 -u 1.25e-8 -g1 -x1 -X1000000 -L -M genome1,genome2, $DIR/plot_sim $DIR/ms2g1.psmc $DIR/ms2g2.psmc
./utils/MS2JSFS.py <(gunzip -c $DIR/sim.ms.gz) ms2g1 ms2g2 > $DIR/sim.jsfs
if [ $CLEAN -eq 1 ]; then
#	head -n1 $DIR/sim.ms > $DIR/command.ms
	rm $DIR/sim.ms.gz
	rm $DIR/ms2g1.ms
	rm $DIR/ms2g2.ms
fi
#./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR -tol 1e-4 -sm 70 -o solution.migr > $DIR/output.txt
#./MigrationPlot.py ms2g1.psmc ms2g2.psmc solution.migr -wd $DIR
