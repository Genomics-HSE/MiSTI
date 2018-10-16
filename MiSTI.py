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

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('fpsmc1',
                    help='psmc file 1')
parser.add_argument('fpsmc2',
                    help='psmc file 2')
parser.add_argument('fjafs',
                    help='joint allele frequency spectrum file')

parser.add_argument('-o', '--fout', nargs=1, default='',
                    help='output file, default is stdout')
parser.add_argument('-wd', nargs=1, default='',
                    help='working directory (path to data files)')
parser.add_argument('-pr', nargs=1, type=int, default=1,
                    help='number of processes for multiprocessing optimisation (default is 1)')
parser.add_argument('-tol', nargs=1, type=float, default=1e-4,
                    help='optimisation precision (default is 1e-4)')
                    
parser.add_argument('-mu0', nargs=2, type=float, default=[0.0, 0.0],
                    help='initial values for mu0 in the optimisation')

parser.add_argument('-ms', '--migstart', nargs=1, type=int, default=0,
                    help='time when migration starts (migration is zero in the recent past)')
parser.add_argument('-me', '--migend', nargs=1, type=int, default=0,
                    help='time when migration ends (migration is zero in the distant past)')

parser.add_argument('-sm', nargs=1, type=float, default=0,
                    help='minimal split time')
parser.add_argument('-sM', nargs=1, type=float, default=0,
                    help='maximal split time')
parser.add_argument('-fil', nargs=1, type=int, default=0, #fil for first interval length in case of optimisation of first lambdas
                    help='first interval length')
parser.add_argument('--sdate', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome)') 
parser.add_argument('--discr', '-d', nargs=1, type=int, default=1,
                    help='discritesation of intervals (default is 1 for no discritisation).')
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')
parser.add_argument('-oml', action='store_true',
                    help='Optimisation of migration rates and lambdas')
parser.add_argument('-ol', action='store_true',
                    help='Optimisation of lambdas')

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
if isinstance(clargs.sm, list):
    clargs.sm = clargs.sm[0]
if isinstance(clargs.sM, list):
    clargs.sM = clargs.sM[0]
if isinstance(clargs.sdate, list):
    clargs.sdate = clargs.sdate[0]
if isinstance(clargs.fil, list):
    clargs.fil = clargs.fil[0]
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
if isinstance(clargs.oml, list):
    clargs.oml = clargs.oml[0]
if isinstance(clargs.ol, list):
    clargs.ol = clargs.ol[0]
if isinstance(clargs.uf, list):
    clargs.uf = clargs.uf[0]
if isinstance(clargs.migstart, list):
    clargs.migstart = clargs.migstart[0]
if isinstance(clargs.migend, list):
    clargs.migend = clargs.migend[0]
if isinstance(clargs.funits, list):
    clargs.funits = clargs.funits[0]

if isinstance(clargs.debug, list):
    clargs.debug = clargs.debug[0]

if clargs.migstart >= clargs.migend and clargs.migend > 0:
    PrintErr("--migend should be larger than --migstart. Aborted.")
    sys.exit(0)

if clargs.debug:
    try:
        from guppy import hpy
        loaded_guppy = True
        h = hpy()
    except ImportError:
        loaded_guppy = False
        PrintErr("Cannot import hpy from guppy module. Profiling is disabled.")

mode = "optimize"
if clargs.llh:
    mode = "llhmodel"

clargs.settings = {
    "enableOutput": False,
    "nosmooth": False,
    "unfolded": clargs.uf,
    "sampleDate": clargs.sdate
}

units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)

def Optimize(times, lambdas, dataJAFS):
    global clargs
    PrintErr("Number of processes: ", clargs.pr)
    smin = int(clargs.sm)
    smax = int(clargs.sM)
    PrintErr("Optimizing for split time range from ", smin, " to ", smax)
    PrintErr("Optimization tolerance ", clargs.tol)
    splitTimes = list(range( smin, smax ))
    res = []
    if clargs.oml:
        mu0 = clargs.mu0+[1.0, 1.0]
        data = [times, lambdas, dataJAFS, 0, mu0]
        for splitT in splitTimes:
            data[3] = splitT
            data[4] = mu0
            res.append( RunSolveMuLa(data) )
            mu0 = res[-1][0]
    elif clargs.ol:
        mu0 = clargs.mu0+[1.0, 1.0]
        data = [times, lambdas, dataJAFS, 0, mu0]
        for splitT in splitTimes:
            data[3] = splitT
            data[4] = mu0
            res.append( RunSolveLa(data) )
            mu0 = res[-1][0]
    else:
        mu0 = clargs.mu0
        if clargs.pr == 1:
            data = [times, lambdas, dataJAFS, 0, mu0, ]
            for splitT in splitTimes:
                data[3] = splitT
                data[4] = mu0
                res.append( RunSolve(data) )
                mu0 = res[-1][0]
        else:
            splitVals = [ [times, lambdas, dataJAFS, splitT, mu0] for splitT in splitTimes ]
            p = multiprocessing.Pool( clargs.pr )
            res = p.map(RunSolve, splitVals)
            p.close()
            p.join()
    print(res)
    ct = 0
    for el in res:
        ct += el[3]
    print("Multiprocessing CPU time: ", ct)
    res = sorted( res, key=lambda val: val[1])[-1]
    return(res)

def RunSolve(args):
    global clargs
    t1 = time.process_time()
    PrintErr("Solving for split times ", args[3], ", initial conditions ", args[4])
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], enableOutput = False, smooth = True, trueEPS = clargs.trueEPS, unfolded = clargs.uf, migStart = clargs.migstart, migEnd = clargs.migend)
    muSol = Migration.Solve(clargs.tol, args[4])
    muSol.append(args[3])
    muSol[1] = Migration.JAFSLikelyhood( muSol[0] )
    print( muSol[1] )
    print( muSol )
    MigrationInference.Report()
    t2 = time.process_time()
    muSol.append(t2-t1)
    return( muSol )

def RunSolveMuLa(args):
    global clargs
    t1 = time.process_time()
    PrintErr("Solving for split times ", args[3], ", initial conditions ", args[4])
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = (not clargs.nosmooth), trueEPS = clargs.trueEPS, unfolded = clargs.uf, sampleDate = args[6], migStart = clargs.migstart, migEnd = clargs.migend)
    muSol = Migration.SolveMuLa(clargs.tol, args[4], clargs.fil)
    muSol.append(args[3])
    print(Migration.JAFSLikelyhood( muSol[0] ) )
    MigrationInference.Report()
    t2 = time.process_time()
    muSol.append(t2-t1)
    return( muSol )
    
def RunSolveLa(args):
    global clargs
    t1 = time.process_time()
    PrintErr("Solving for split times ", args[3], ", initial conditions ", args[4])
    Migration = MigrationInference(args[0], args[1], args[2], args[4][0:2], args[3], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = (not clargs.nosmooth), trueEPS = clargs.trueEPS, unfolded = clargs.uf, sampleDate = args[6], migStart = clargs.migstart, migEnd = clargs.migend)
    muSol = Migration.SolveLa(clargs.tol, args[4][2:], clargs.fil)
    muSol.append(args[3])
    print(Migration.JAFSLikelyhood( muSol[0] ) )
    MigrationInference.Report()
    t2 = time.process_time()
    muSol.append(t2-t1)
    return( muSol )

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
smin = min( clargs.sm, len(times) )
smin = max(smin, inputData[6])
smax = min( clargs.sM, len(times)+1 )
smax = max( smax, smin )
if clargs.sM == 0:
    smax = len(times)
if smax == smin:
    print("-sm should be strictly smaller than -sM.")
    sys.exit(0)
clargs.sm = smin
clargs.sM = smax
print("Split time will be varied from", smin, "up to", smax)


sol = [[], [], []]
if mode == "optimize":
    sol = Optimize(inputData[0], inputData[1], dataJAFS)
    print(sol)
elif mode == "llhmodel":
    res = []
    for splitT in range( math.floor(clargs.sm), math.ceil(clargs.sM) ):
        discr = clargs.discr
        for ds in range(discr):
            sT = splitT + ds/discr
            if sT < clargs.sm or sT >= clargs.sM:
                continue
            Migration = MigrationInference(inputData[0][:], inputData[1][:], dataJAFS, clargs.mu0, sT, thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = True, unfolded = clargs.uf, trueEPS = clargs.trueEPS, migStart = clargs.migstart, migEnd = clargs.migend, sampleDate = inputData[6])
            llh_tmp = Migration.JAFSLikelyhood( clargs.mu0 )
            print("splitT = ", sT, "\tlikelihood = ", llh_tmp)
            res.append( [clargs.mu0, llh_tmp, sT, 0.0] )
    if clargs.debug:
        print(res)
    sol = sorted( res, key=lambda val: val[1])[-1]
    confInt = [None, None]
    confIntTmp = [val for val in res if val[1] <= sol[1]-1.92 and val[2] < sol[2]]
    if len(confIntTmp) == 0:
        confInt[0] = sol
    else:
        confInt[0] = sorted( confIntTmp, key=lambda val: val[2])[-1]
    confIntTmp = [val for val in res if val[1] <= sol[1]-1.92 and val[2] > sol[2]]
    if len(confIntTmp) == 0:
        confInt[1] = sol
    else:
        confInt[1] = sorted( confIntTmp, key=lambda val: val[2])[0]
    if clargs.debug:
        print(confInt)
    
else:
    PrintErr("Unknown mode")
    sys.exit(0)

print("\n\nParameter estimates:")

splitT = sol[2]
print("splitT = ", splitT, "\ttime = ", (sum(inputData[0][0:int(splitT)])+inputData[0][int(splitT)]*(splitT%1))*inputData[2], "\tmigration rates = ", sol[0][0]/migrUnit, ", ", sol[0][1]/migrUnit, "\tllh = ", sol[1])
#print("\tmigStart = ", clargs.migstart, "\tmigration start time = ", sum(inputData[0][0:clargs.migstart])*inputData[2])
if mode == "llhmodel":
    print("Confidence interval: ", confInt[0][2] , " ", confInt[1][2], "\t", (sum(inputData[0][0:int(confInt[0][2])])+inputData[0][int(confInt[0][2])]*(confInt[0][2]%1))*inputData[2], "\t", (sum(inputData[0][0:int(confInt[1][2])])+inputData[0][int(confInt[1][2])]*(confInt[1][2]%1))*inputData[2])

print("\n")
Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, sol[0], sol[2], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = True, unfolded = clargs.uf, trueEPS = clargs.trueEPS, migStart = clargs.migstart, migEnd = clargs.migend, sampleDate = inputData[6])
migrationIO.OutputMigration(fout, sol[0], Migration)

#MigrationInference.Report()
t2 = time.time()
if clargs.debug:
    PrintErr("Total time ", t2-t1)
print("Total time ", t2-t1)
if clargs.debug:
    if loaded_guppy:
        print(h.heap())
sys.exit(1)
