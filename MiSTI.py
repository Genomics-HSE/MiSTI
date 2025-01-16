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
import random
from math import (exp,log, ceil)
import time
import multiprocessing
from MigrationInference import MigrationInference
from CorrectLambda import CorrectLambda
import migrationIO
from migrationIO import PrintErr

random.seed()
t0 = time.time()

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('fpsmc1',
                    help='psmc file 1')
parser.add_argument('fpsmc2',
                    help='psmc file 2')
parser.add_argument('fjafs',
                    help='joint allele frequency spectrum file')
parser.add_argument('st', nargs=1, type=float,
                    help='split time')

parser.add_argument('-o', '--fout', nargs=1, default='',
                    help='output file, default is stdout')
parser.add_argument('-wd', nargs=1, default='',
                    help='working directory (path to data files)')
parser.add_argument('-tol', nargs=1, type=float, default=1e-4,
                    help='optimisation precision (default is 1e-4)')
parser.add_argument('-mth', nargs=1, type=float, default=0.0,
                    help='mixture treshhold (default is 0.0)')

parser.add_argument('-mi', nargs=5, action = 'append',
                    help='migration rate, require 5 arguments:\n\t\tsource population index (1 or 2)\n\t\tmigration start time\n\t\tmigration end time\n\t\tmigration rate initial value\n\t\tfixed(0) or optimised(1) parameter.')#-mi [npop:1/2] [migStart] [migEnd] [init val] [var:0/1]
parser.add_argument('-pu', nargs=4, action = 'append',
                    help='pulse migration, require 4 arguments:\n\t\tsource population index (1 or 2)\n\t\tpulse migration time\n\t\tpulse migration rate\n\t\tfixed(0) or optimised(1) parameter.')#-mi [npop:1/2] [migStart] [migEnd] [init val] [var:0/1]

parser.add_argument('--sdate', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome)')
parser.add_argument('--hetloss', '-hl', nargs=2, type=float,
                    help='loss of heterozygosity for the first and the second genomes (default is 0.0)')
parser.add_argument('--discr', '-d', nargs=1, type=int, default=1,
                    help='discritesation of intervals (default is 1 for no discritisation).')
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')


parser.add_argument('-uf', action='store_true',
                    help='Unfolded spectrum')
parser.add_argument('--nosmooth', action='store_true',
                    help='Don\'t smooth (don\'t make constant on the psmc time intervals)')
parser.add_argument('--trueEPS', action='store_true',
                    help='Consider input as true effective population size (instead of mixed coalescence rates)')
parser.add_argument('--cpfit', action='store_true',
                    help='Approximate EPS by fitting probabilities to coalesce within each interval (default is fitting expected coalescence time within each interval)')
#parser.add_argument('--bsSize', '-bs', type=int, default=0,
#                    help='Number of bootstrap repetitions')
parser.add_argument('--bsMode', '-bs', nargs=1, type=int, default=-1,
                    help='Generate single bootstrap sample')


parser.add_argument('--psmcMode', '-pm', type=int, default=0,
                    help='PSMC mode')

parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')


clargs = parser.parse_args()

if isinstance(clargs.fout, list):
    clargs.fout = clargs.fout[0]
if isinstance(clargs.wd, list):
    clargs.wd = clargs.wd[0]
if isinstance(clargs.tol, list):
    clargs.tol = clargs.tol[0]
if isinstance(clargs.mth, list):
    clargs.mth = clargs.mth[0]
if isinstance(clargs.st, list):
    clargs.st = clargs.st[0]
if isinstance(clargs.sdate, list):
    clargs.sdate = clargs.sdate[0]
if isinstance(clargs.discr, list):
    clargs.discr = clargs.discr[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]
if isinstance(clargs.nosmooth, list):
    clargs.nosmooth = clargs.nosmooth[0]
if isinstance(clargs.trueEPS, list):
    clargs.trueEPS = clargs.trueEPS[0]
if isinstance(clargs.cpfit, list):
    clargs.cpfit = clargs.cpfit[0]
if isinstance(clargs.uf, list):
    clargs.uf = clargs.uf[0]
if isinstance(clargs.bsMode, list):
    clargs.bsMode = clargs.bsMode[0]
if isinstance(clargs.psmcMode, list):
    clargs.psmcMode = clargs.psmcMode[0]
if clargs.mi is None:
    clargs.mi=[]
if clargs.pu is None:
    clargs.pu=[]

if isinstance(clargs.funits, list):
    clargs.funits = clargs.funits[0]

if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

clargs.settings = {
    "enableOutput": False,
    "unfolded": clargs.uf,
    "sampleDate": clargs.sdate
}

units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)
units.PrintUnits()
if clargs.hetloss is not None:
    units.SetHetLoss(clargs.hetloss)

print( " ".join(sys.argv) )

startTime = time.strftime("Job run at %H:%M:%S on %d %b %Y")
if clargs.debug:
    PrintErr(startTime)
print(startTime)

fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc1 )
fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc2 )
fjafs  = os.path.join( clargs.wd, clargs.fjafs  )

if clargs.debug:
    PrintErr("Reading from files: ", fpsmc1, ", ", fpsmc2, ", ", fjafs)
print("Reading from files:")
print("pop1\t", fpsmc1)
print("pop2\t", fpsmc2)
print("jafs\t", fjafs)

dataJAFS = migrationIO.ReadJAFS(fjafs)
if clargs.bsMode == -1:
    inputSFS = [0 for _ in range(8)]
    for sfs in dataJAFS.jafs:
        inputSFS = [v+u for v, u in zip(inputSFS, sfs)]
else:
    inputSFS = dataJAFS.jafs[clargs.bsMode]


print("IMPORTANT NOTICE!!! Every time you are running MiSTI, make sure that psmc files are supplied in the same order as populations appear in the joint allele frequency spectrum.")

fout   = clargs.fout
if fout != "":
    fout  = os.path.join( clargs.wd, clargs.fout  )

if clargs.debug:
    print(clargs)

if clargs.psmcMode == 0:
    inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2, clargs.sdate, clargs.rd)
else:
    inputData = migrationIO.ReadPSMC1(fpsmc1, fpsmc2, clargs.rd, divergenceTime = clargs.st)
if inputData.divergenceTime == -1:
    inputData.divergenceTime = clargs.st

if clargs.debug:
    print("INPUT DATA")
    for v in inputData:
        print(v)

    print("END INPUT DATA")

migrUnit = inputData.scaleEPS/2#Convert to ms migration rates (1/2 factor!)

times = inputData.times


sol = [[], [], []]

t1 = time.time()

Migration = MigrationInference(inputData.times, inputData.lambdas, inputSFS, inputData.divergenceTime, clargs.mi, clargs.pu, thrh = [inputData.theta, inputData.rho], Tpsmc = inputData.Tpsmc, enableOutput = False, smooth = not clargs.nosmooth, unfolded = clargs.uf, trueEPS = clargs.trueEPS, sampleDate = inputData.sampleDateDiscr, mixtureTH = clargs.mth, cpfit = clargs.cpfit)
sol = Migration.Solve(clargs.tol)
print(sol)


print("\nParameter estimates:")

migStr = "migration rates "
migFixed = []
for el in clargs.mi:
    if int(el[4]) == 0:
        migFixed.append(float(el[3]))
if len(migFixed) > 0:
    migFixedStr = "fixed = [" + ", ".join([str(v) for v in migFixed]) + "]"
else:
    migFixedStr = ""

if len(sol[0]) > 0:
    migOptStr = "optim = [" + ", ".join([str(v) for v in sol[0]]) + "]"
else:
    migOptStr = ""

if migFixedStr != "" and migOptStr != "":
    migStr = migFixedStr + "\t" + migOptStr
else:
    migStr = migFixedStr + migOptStr
#print("splitT =", inputData.divergenceTime, "\ttime =", (sum(inputData[0][0:int(splitT)])+inputData[0][int(splitT)]*(splitT%1))*inputData[2], "\tmigration rates =", [v/migrUnit for v in sol[0]], "\tllh =", sol[1])
print("bs_id =", clargs.bsMode, "\tsplitT =", inputData.divergenceTime, "\ttime =", sum(inputData.times[0:ceil(inputData.divergenceTime)])*inputData.scaleTime, "\tmigration rates", migStr, "\tllh =", sol[1])
print("\n")

t2 = time.time()

if sol[1] == -10**9:
    print("Failed to fit such a model.")
else:
    if clargs.bsMode == 0:
        migrationIO.OutputMigration(fout, sol[0], Migration, inputData.scaleTime, inputData.scaleEPS)

t3 = time.time()

MigrationInference.Report()
if clargs.debug:
    PrintErr("Runtime:   optimisation ", t2-t1)
    PrintErr("           total        ", t3-t0)

print("Runtime:   optimisation", t2-t1)
print("           total       ", t3-t0)
sys.exit(0)
