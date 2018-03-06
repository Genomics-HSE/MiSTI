#!/usr/bin/env python3

import sys
import collections
import numpy
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)
import time
import matplotlib.pyplot as plt
import multiprocessing
from MigrationInference import MigrationInference
import migrationIO


def Optimize(times, lambdas, dataJAFS, procNum=4):
    p = multiprocessing.Pool(procNum)
    splitVals = range( len(times) )
    splitVals = [ [times, lambdas, dataJAFS, splitT] for splitT in [98,99,100,101] ]
    res = p.map(RunSolve, splitVals)
    print(res)

#inputData[0], inputData[1], dataJAFS, splitT
def RunSolve(args):
    print("Solving for split times ", args[3])
    Migration = MigrationInference(args[0], args[1], args[2], [0,0], args[3], 1.0, correct = True, enableOutput = False, smooth = True)
    muSol = Migration.Solve()
    print(muSol)


if len(sys.argv) < 4:
    print("./migration <PSMC input file 1> <PSMC input file 2> <JAF spectrum file>")
    sys.exit(0)
t1 = time.clock()
fpsmc1 = sys.argv[1]
fpsmc2 = sys.argv[2]
fjafs  = sys.argv[3]
doPlot = True
inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2, doPlot = doPlot, skip = 0)
dataJAFS = migrationIO.ReadJAFS(fjafs)




maxllh = 0.0
maxmu = []
maxsplitT = 0

Optimize(inputData[0], inputData[1], dataJAFS)



if doPlot:
    mu = maxmu
    splitT = maxsplitT
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, correct = True, smooth = True)
    Migration.JAFSLikelyhood( mu )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times))]
    print( Migration.lc[0:-1] )
    print(mu)
    plt.step([v*inputData[2] for v in times[0:]], [1.0/max(v[0],0.1)*inputData[3] for v in Migration.lc[0:-1]])
    splT=sum(inputData[0][0:splitT])
    plt.axvline(splT*inputData[2], color='r')
    
if doPlot:
    plt.savefig("data_bn20_bn20/plot.png")

MigrationInference.Report()

t2 = time.clock()
print("Total time ", t2-t1)


#DIR=
#./migration.py $DIR/ms2g1.psmc $DIR/ms2g2.psmc $DIR/sim.jafs >> $DIR/output.txt