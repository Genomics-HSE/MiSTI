#!/usr/bin/env python3

#    Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)
#
#    This file is part of MiSTI.
#
#    MiSTI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MiSTI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MiSTI.  If not, see <https://www.gnu.org/licenses/>.

import sys
import os

os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import collections
import argparse
import numpy
import math
from math import (exp,log)
import time
import multiprocessing
from MigrationInference import MigrationInference
from CorrectLambda import CorrectLambda
import migrationIO
from migrationIO import PrintErr

t0 = time.time()

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('msstring',
                    help='ms style command')
parser.add_argument('fjafs', nargs='?', default='',
                    help='joint allele frequency spectrum file')

parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')

parser.add_argument('-uf', action='store_true',
                    help='Unfolded spectrum')
parser.add_argument('--bsSize', '-bs', type=int, default=0,
                    help='Number of bootstrap repetitions')

parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')


clargs = parser.parse_args()

if isinstance(clargs.uf, list):
    clargs.uf = clargs.uf[0]
if isinstance(clargs.bsSize, list):
    clargs.bsSize = clargs.bsSize[0]

units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)
units.PrintUnits()

jafs_input = False
if clargs.fjafs == "":
    inputSFS = [0 for _ in range(8)]
else:
    dataJAFS = migrationIO.ReadJAFS(clargs.fjafs)
    jafs_input = True
    inputSFS = [0 for _ in range(8)]
    for sfs in dataJAFS.jafs:
        inputSFS = [v+u for v, u in zip(inputSFS, sfs)]

inputData = migrationIO.ReadMS(clargs.msstring)
if False:
    print("INPUT DATA")
    for v in inputData:
        print(v)
    print("END INPUT DATA")

Migration = MigrationInference(inputData[0], inputData[1], inputSFS, inputData[2], inputData[3], inputData[4], unfolded = clargs.uf, trueEPS = True)
llh = Migration.JAFSLikelihood([])
print("Expected SFS", Migration.JAFS)
if jafs_input:
    jafs = inputSFS[1:]
    norm = sum(jafs)
    jafs = [v/norm for v in jafs]
    print("Data     SFS", jafs)
    print("data llh under the model is", llh)
    mllh = Migration.MaximumLLHFunction()
    print("maximum of the llh function is", mllh)

    if clargs.bsSize > 1:
        bs_llh = []
        bs_size = clargs.bsSize
        for i in range(clargs.bsSize):
            Migration.SetJAFS(migrationIO.BootstrapJAFS(dataJAFS))
            bs_llh.append( Migration.JAFSLikelihood(sol[0]) )
        bs_llh.sort()
        cutoff = math.ceil(0.05*bs_size)
        print("10% confidence interval", bs_llh[cutoff], bs_llh[-cutoff])
        cutoff = math.ceil(0.025*bs_size)
        print("5% confidence interval", bs_llh[cutoff], bs_llh[-cutoff])


intervals = []
ct = 0
for i in range(Migration.numT):
    intervals.append([ct+Migration.times[i], Migration.lh[i][0], Migration.lh[i][1], , Migration.mi[i][0], Migration.mi[i][1]])
Migration.CoalescentRates()#interval = [time, lambda1, lambda2, mu1, mu2], discr = number of intervals in the discretization
migrationIO.OutputMigration(fout, sol[0], Migration)
    
#print("splitT =", splitT, "\ttime =", (sum(inputData[0][0:int(splitT)])+inputData[0][int(splitT)]*(splitT%1))*inputData[2], "\tmigration rates =", [v/migrUnit for v in sol[0]], "\tllh =", sol[1])


sys.exit(1)
