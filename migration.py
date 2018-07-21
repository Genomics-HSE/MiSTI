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

parser.add_argument('-sm', nargs=1, type=int, default=0,
                    help='minimal split time')
parser.add_argument('-sM', nargs=1, type=int, default=0,
                    help='maximal split time')
parser.add_argument('-fil', nargs=1, type=int, default=0, #fil for first interval length in case of optimisation of first lambdas
                    help='first interval length')
parser.add_argument('-sd', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome)')
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('-oml', action='store_true',
                    help='Optimisation of migration rates and lambdas')
parser.add_argument('-ol', action='store_true',
                    help='Optimisation of lambdas')

parser.add_argument('-uf', action='store_true',
                    help='Unfolded spectrum')
parser.add_argument('-llh', action='store_true',
                    help='Compute model llh')
parser.add_argument('--smooth', action='store_true',
                    help='Smooth (make constant on the psmc time intervals)')
parser.add_argument('--trueEPS', action='store_true',
                    help='Consider input as true effective population size (instead of mixed coalescence rates)')


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
if isinstance(clargs.sM, list):
    clargs.sd = clargs.sd[0]
if isinstance(clargs.fil, list):
    clargs.fil = clargs.fil[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]
if isinstance(clargs.llh, list):
    clargs.llh = clargs.llh[0]
if isinstance(clargs.smooth, list):
    clargs.smooth = clargs.smooth[0]
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

mode = "optimize"
if clargs.llh:
    mode = "llhmodel"

clargs.settings = {
    "enableOutput": False,
    "smooth": True,
    "unfolded": clargs.uf,
    "sampleDate": clargs.sd
}

def Optimize(times, lambdas, dataJAFS):
    global clargs
    PrintErr("Number of processes: ", clargs.pr)
    smin = clargs.sm
    smax = clargs.sM
    '''smin = min( clargs.sm, len(times) )
    smax = min( clargs.sM, len(times)+1 )
    smax = max( smax, smin )
    print("smax = ", smax, "\tsmin = ", smin)
    if clargs.sM == 0:
        smax = len(times)
    if smax == smin:
        print("-sm should be strictly smaller than -sM.")
        sys.exit(0)'''
    PrintErr("Optimizing for split time range from ", smin, " to ", smax)
    PrintErr("Optimization tolerance ", clargs.tol)
    splitTimes = list(range( smin, smax ))
#    splitTimes = list( range(100, 102) )
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
            data = [times, lambdas, dataJAFS, 0, mu0]
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
        
#    p.close()
#    res = sorted( res, key=lambda val: val[2])
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
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = clargs.smooth, trueEPS = clargs.trueEPS, unfolded = clargs.uf, sampleDate = clargs.sd, migStart = clargs.migstart, migEnd = clargs.migend)
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
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = clargs.smooth, trueEPS = clargs.trueEPS, unfolded = clargs.uf, sampleDate = clargs.sd, migStart = clargs.migstart, migEnd = clargs.migend)
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
    Migration = MigrationInference(args[0], args[1], args[2], args[4][0:2], args[3], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = clargs.smooth, trueEPS = clargs.trueEPS, unfolded = clargs.uf, sampleDate = clargs.sd, migStart = clargs.migstart, migEnd = clargs.migend)
    muSol = Migration.SolveLa(clargs.tol, args[4][2:], clargs.fil)
    muSol.append(args[3])
    print(Migration.JAFSLikelyhood( muSol[0] ) )
    MigrationInference.Report()
    t2 = time.process_time()
    muSol.append(t2-t1)
    return( muSol )


print( " ".join(sys.argv) )
#t1 = time.clock()
t1 = time.time()

startTime = time.strftime("Job run at %H:%M:%S on %d %b %Y")
PrintErr(startTime)
print(startTime)

fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc1 )
fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc2 )
fjafs  = os.path.join( clargs.wd, clargs.fjafs  )

'''cl = CorrectLambda()
mu1=[1.0,1.0]
inter1 = [0, 1, 1, mu1[0], mu1[1]]
inter2 = [0.02, 0.05, 0.05, mu1[0], mu1[1]]
inter3 = [0.075, 0.5, 0.5, mu1[0], mu1[1]]
inter4 = [2.5, 1, 1, mu1[0], mu1[1]]
splitT = 2.5
intervals = [inter1, inter2, inter3, inter4]
coalRates = cl.CoalRates(intervals, splitT, 100)
print(coalRates)
sys.exit(0)'''


PrintErr("Reading from files: ", fpsmc1, ", ", fpsmc2, ", ", fjafs)
print("Reading from files: ", fpsmc1, ", ", fpsmc2, ", ", fjafs)

fout   = clargs.fout
if fout != "":
    fout  = os.path.join( clargs.wd, clargs.fout  )
    
print(clargs)

inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2, clargs.rd)

if 1:
    timesMS = [0, 0.0275, 0.0475, 0.175, 0.75, 3.75, 10]
    epsMS = [[13, 0.25], [0.5, 0.4], [0.5, 0.5], [3, 3], [2, 2], [3, 3], [6, 6]]
    inputData[0] = [2*(u-v) for u, v in zip(timesMS[1:], timesMS[:-1])]
    inputData[1] = [[1/u[0], 1/u[1]] for u in epsMS]
#-n 1 13 -n 2 0.25
#-en 0.0275 1 0.5 -em 0.0275 2 1 12.5 -ej 0.0475 2 1 -eM 0.0475 0.0
#-eN 0.175 3
#-eN 0.75 2
#-eN 3.75 3
#-eN 10 6
#print(inputData)
migrUnit = inputData[3]/2#Convert to ms migration rates (1/2 factor!)
#migrUnit = (inputData[1][0][0]+inputData[1][1][0])/4.0
dataJAFS = migrationIO.ReadJAFS(fjafs)

times = inputData[0]
smin = min( clargs.sm, len(times) )
smax = min( clargs.sM, len(times)+1 )
smax = max( smax, smin )
if clargs.sM == 0:
    smax = len(times)
if smax == smin:
    print("-sm should be strictly smaller than -sM.")
    sys.exit(0)
clargs.sm = smin
clargs.sM = smax


sol = [[], [], []]
if mode == "optimize":
    sol = Optimize(inputData[0], inputData[1], dataJAFS)
    print(sol)
elif mode == "llhmodel":
    res = []
    for splitT in range( clargs.sm, clargs.sM ):
        Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, clargs.mu0, splitT, thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = clargs.smooth, unfolded = clargs.uf, trueEPS = clargs.trueEPS, migStart = clargs.migstart, migEnd = clargs.migend)
        llh_tmp = Migration.JAFSLikelyhood( clargs.mu0 )
        print("splitT = ", splitT, "\tlikelihood = ", llh_tmp)
        res.append( [clargs.mu0, llh_tmp, splitT, 0.0] )
    print(res)
    sol = sorted( res, key=lambda val: val[1])[-1]
else:
    PrintErr("Unknown mode")
    sys.exit(0)
#sol[0] = [0.34646987, 0.32497276]
#sol[2] = 110

splitT = sol[2]
print("splitT = ", splitT, "\ttime = ", sum(inputData[0][0:splitT])*inputData[2], "\tmu = ", [sol[0][0]/migrUnit,sol[0][1]/migrUnit], "\tllh = ", sol[1])
print("\tmigStart = ", clargs.migstart, "\tmigration start time = ", sum(inputData[0][0:clargs.migstart])*inputData[2])
print("N_0 = ")
Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, sol[0], sol[2], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = clargs.smooth, unfolded = clargs.uf, trueEPS = clargs.trueEPS, migStart = clargs.migstart, migEnd = clargs.migend)
migrationIO.OutputMigration(fout, sol[0], Migration)

#MigrationInference.Report()
#t2 = time.clock()
t2 = time.time()
PrintErr("Total time ", t2-t1)
print("Total time ", t2-t1)
sys.exit(0)

inputData = migrationIO.ReadMS("4 1000 -t 8196 -r 1355 3000000 -l -I 2 2 2 -n 2 1.0 -em 0.0 1 2 2.0 -em 0.0 2 1 2.0 -en 0.01 1 0.05 -en 0.01 2 0.05 -en 0.0375 1 0.5 -en 0.0375 2 0.5 -ej 1.25 2 1 -eM 1.25 0.0 -eN 1.25 1.0")
Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, [1, 1], inputData[4], thrh = [inputData[4], inputData[5]], enableOutput = False, smooth = False, unfolded = clargs.uf, correct = False, migStart = clargs.migstart, migEnd = clargs.migend)
print(Migration.JAFSLikelyhood( [1.0, 1.0] ) )
print(Migration.JAFSLikelyhood( [1.0036919845350205, 1.0035904582181976] ) )
sys.exit(0)

#DIR=
#./migration.py $DIR/ms2g1.psmc $DIR/ms2g2.psmc $DIR/sim.jafs >> $DIR/output.txt
#./migration.py ms2g1.psmc ms2g2.psmc sim.jafs -wd $DIR >> $DIR/output1.txt