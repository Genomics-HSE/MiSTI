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
from MigrationInference2 import MigrationInference
from CorrectLambda import CorrectLambda
import migrationIO
from migrationIO import PrintErr

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
parser.add_argument('-pr', nargs=1, type=int, default=1,
                    help='number of processes for multiprocessing optimisation (default is 1)')
parser.add_argument('-tol', nargs=1, type=float, default=1e-4,
                    help='optimisation precision (default is 1e-4)')
                    
parser.add_argument('-mi', nargs=5, action = 'append',
                    help='initial values for mu0 in the optimisation')#-mi [npop:1/2] [migStart] [migEnd] [init val] [var:0/1]

parser.add_argument('--sdate', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome)') 
parser.add_argument('--discr', '-d', nargs=1, type=int, default=1,
                    help='discritesation of intervals (default is 1 for no discritisation).')
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')


parser.add_argument('-uf', action='store_true',
                    help='Unfolded spectrum')
parser.add_argument('-llh', action='store_true',
                    help='Compute model llh')
parser.add_argument('--nosmooth', action='store_false',
                    help='Don\'t smooth (make constant on the psmc time intervals)')
parser.add_argument('--trueEPS', action='store_true',
                    help='Consider input as true effective population size (instead of mixed coalescence rates)')


parser.add_argument('--debug', action='store_true',
                    help='Debug mode, more input enabled')


clargs = parser.parse_args()

if isinstance(clargs.fout, list):
    clargs.fout = clargs.fout[0]
if isinstance(clargs.wd, list):
    clargs.wd = clargs.wd[0]
if isinstance(clargs.pr, list):
    clargs.pr = clargs.pr[0]
if isinstance(clargs.tol, list):
    clargs.tol = clargs.tol[0]
if isinstance(clargs.st, list):
    clargs.st = clargs.st[0]
if isinstance(clargs.sdate, list):
    clargs.sdate = clargs.sdate[0]
if isinstance(clargs.discr, list):
    clargs.discr = clargs.discr[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]
if isinstance(clargs.llh, list):
    clargs.llh = clargs.llh[0]
if isinstance(clargs.nosmooth, list):
    clargs.nosmooth = clargs.nosmooth[0]
if isinstance(clargs.trueEPS, list):
    clargs.trueEPS = clargs.trueEPS[0]
if isinstance(clargs.uf, list):
    clargs.uf = clargs.uf[0]
if clargs.mi == None:
    clargs.mi=[]

if isinstance(clargs.funits, list):
    clargs.funits = clargs.funits[0]

if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

mode = "optimize"
if clargs.llh:
    mode = "llhmodel"

clargs.settings = {
    "enableOutput": False,
    "nosmooth": False,
    "unfolded": clargs.uf,
    "sampleDate": clargs.sdate
}

print(clargs)

units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)

print( " ".join(sys.argv) )
t1 = time.time()

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
print("IMPORTANT NOTICE!!! Every time you are running MiSTI, make sure that psmc file are supplied in the same order as populations appear in the joint allele frequency spectrum.")

fout   = clargs.fout
if fout != "":
    fout  = os.path.join( clargs.wd, clargs.fout  )

if clargs.debug:
    print(clargs)

inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2, clargs.sdate, clargs.rd)
if clargs.debug:
    print("INPUT DATA")
    for v in inputData:
        print(v)
        
    print("END INPUT DATA")

migrUnit = inputData[3]/2#Convert to ms migration rates (1/2 factor!)

times = inputData[0]


sol = [[], [], []]
if mode == "optimize":
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, clargs.mi, clargs.st, thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = True, unfolded = clargs.uf, trueEPS = clargs.trueEPS, sampleDate = inputData[6])
    sol = Migration.Solve(clargs.tol)
    print(sol)
elif mode == "llhmodel":
    print("Not in this version...")
    sys.exit(0)
    
else:
    PrintErr("Unknown mode")
    sys.exit(0)

print("\n\nParameter estimates:")

splitT = clargs.st

print("splitT = ", splitT, "\ttime = ", (sum(inputData[0][0:int(splitT)])+inputData[0][int(splitT)]*(splitT%1))*inputData[2], "\tmigration rates = ", [v/migrUnit for v in sol[0]], "\tllh = ", sol[1])
print("\n")

migrationIO.OutputMigration2(fout, sol[0], Migration)

#MigrationInference.Report()
t2 = time.time()
if clargs.debug:
    PrintErr("Total time ", t2-t1)
print("Total time ", t2-t1)
sys.exit(1)
