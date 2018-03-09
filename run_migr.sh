#!/bin/sh
if [ "$#" -lt 1 ]; then
  echo "Usage: $0 DIRECTORY [-c]" >&2
  exit 1
fi
DIR=$1

CLEAN=0
if [ "$#" -gt 1 ]; then
	if [ "$2" = "-c" ]; then
		CLEAN=1
	fi
fi

if [ $CLEAN -eq 1 ]; then
	rm $DIR/sim.ms
	rm $DIR/ms2g1.ms
	rm $DIR/ms2g2.ms
fi
./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR -tol 1e-4 -sm 80 -sM 110 -pr 6 -o solution.migr > $DIR/output.txt
DIR=data_bn00_bn05
./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR -tol 1e-4 -sm 80 -sM 110 -pr 6 -o solution.migr > $DIR/output.txt
DIR=data_bn00_bn10
./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR -tol 1e-4 -sm 80 -sM 110 -pr 6 -o solution.migr > $DIR/output.txt
DIR=data_bn00_bn20
./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR -tol 1e-4 -sm 80 -sM 110 -pr 6 -o solution.migr > $DIR/output.txt
#./MigrationPlot.py ms2g1.psmc ms2g2.psmc solution.migr -wd $DIR
