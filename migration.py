#!/usr/bin/env python3

import sys
import os
import collections
import argparse
import numpy
import math
from math import (exp,log)
import time
import multiprocessing
from MigrationInference import MigrationInference
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
parser.add_argument('-tol', nargs=1, type=float, default=1,
                    help='optimisation precision (default is 1e-4)')
parser.add_argument('-sm', nargs=1, type=int, default=0,
                    help='minimal split time')
parser.add_argument('-sM', nargs=1, type=int, default=0,
                    help='maximal split time')

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

def Optimize(times, lambdas, dataJAFS):
    global clargs
    PrintErr("Number of processes: ", clargs.pr)
    smin = min( clargs.sm, len(times) )
    smax = min( clargs.sM, len(times) )
    smax = max( smax, smin )
    if clargs.sM == 0:
        smax = len(times)
    if smax == smin:
        print("-sm should be strictly smaller than -sM.")
        sys.exit(0)
    PrintErr("Optimizing for split time range from ", smin, " to ", smax)
    PrintErr("Optimization tolerance ", clargs.tol)
    splitTimes = list(range( smin, smax ))
#    splitTimes = list( range(100, 102) )
    splitVals = [ [times, lambdas, dataJAFS, splitT] for splitT in splitTimes ]
    res = []
    if clargs.pr == 1:
        for el in splitVals:
            print(el)
            res.append( RunSolve(el) )
    else:
        p = multiprocessing.Pool( clargs.pr )
        res = p.map(RunSolve, splitVals)
        p.close()
        p.join()
#    p.close()
#    res = sorted( res, key=lambda val: val[2])
    print(res)
    res = sorted( res, key=lambda val: val[1])[-1]
    return(res)

def RunSolve(args):
    global clargs
    PrintErr("Solving for split times ", args[3])
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], 1.0, enableOutput = False, smooth = True)
    muSol = Migration.Solve(clargs.tol)
    muSol.append(args[3])
    MigrationInference.Report()
    return( muSol )

#t1 = time.clock()
t1 = time.time()

fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc1 )
fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc2 )
fjafs  = os.path.join( clargs.wd, clargs.fjafs  )

PrintErr("Reading from files: ", fpsmc1, ", ", fpsmc2, ", ", fjafs)
print("Reading from files: ", fpsmc1, ", ", fpsmc2, ", ", fjafs)

fout   = clargs.fout
if fout != "":
    fout  = os.path.join( clargs.wd, clargs.fout  )
    
print(clargs)

inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2)
migrUnit = inputData[1][0][0]/2#TODO remove or fix
dataJAFS = migrationIO.ReadJAFS(fjafs)
sol = [[], [], []]
sol = Optimize(inputData[0], inputData[1], dataJAFS)
#print(sol)
#sol[0] = [0.34646987, 0.32497276]
#sol[2] = 98

splitT = sol[2]
print("splitT = ", splitT, "\ttime = ", sum(inputData[0][0:splitT])*inputData[2], "\tmu = ", [sol[0][0]/migrUnit,sol[0][1]/migrUnit], "\tllh = ", sol[1])
Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, sol[0], sol[2], 1.0, enableOutput = False, smooth = True, correct = True)
migrationIO.OutputMigration(fout, sol[0], Migration)

#MigrationInference.Report()

#t2 = time.clock()
t2 = time.time()
PrintErr("Total time ", t2-t1)
print("Total time ", t2-t1)
sys.exit(0)


#DIR=
#./migration.py $DIR/ms2g1.psmc $DIR/ms2g2.psmc $DIR/sim.jafs >> $DIR/output.txt
#./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR >> $DIR/output1.txt