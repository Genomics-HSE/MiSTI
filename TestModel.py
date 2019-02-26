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

parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')


clargs = parser.parse_args()

#if isinstance(clargs.msstring, list):
#    clargs.msstring = clargs.msstring[0]

units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)
jafs_input = False
if clargs.fjafs == "":
    dataJAFS = migrationIO.JAFS(jafs=[0 for _ in range(8)])
else:
    dataJAFS = migrationIO.ReadJAFS(clargs.fjafs)
    jafs_input = True

inputData = migrationIO.ReadMS(clargs.msstring)
if False:
    print("INPUT DATA")
    for v in inputData:
        print(v)
    print("END INPUT DATA")

Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, inputData[2], inputData[3], thrh = [1, 1], unfolded = clargs.uf, trueEPS = True, mixtureTH = 0)
llh = Migration.JAFSLikelyhood([])
print("Expected SFS", Migration.JAFS)
if jafs_input:
    jafs = dataJAFS.jafs[1:]
    norm = sum(jafs)
    jafs = [v/norm for v in jafs]
    print("Data     SFS", jafs)
    print("data llh under the model is", llh)
    mllh = Migration.MaximumLLHFunction()
    print("maximum of the llh function is", mllh)


#print("splitT =", splitT, "\ttime =", (sum(inputData[0][0:int(splitT)])+inputData[0][int(splitT)]*(splitT%1))*inputData[2], "\tmigration rates =", [v/migrUnit for v in sol[0]], "\tllh =", sol[1])


sys.exit(1)
