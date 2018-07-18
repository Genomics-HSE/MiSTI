#!/usr/bin/env python3

import sys
import os
import collections
import argparse
import numpy
import math
from math import (exp,log)
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import migrationIO
plt.ioff()

parser = argparse.ArgumentParser(description='Migration inference from PSMC help scripts.')
parser.add_argument('fmigr',
                    help='migr file')
clargs = parser.parse_args()
#fmigr  = os.path.join( clargs.wd, clargs.fmigr  )

data = migrationIO.ReadMigration(clargs.fmigr, False)
#data = MigData(splitT = data[0], migStart = data[1], migEnd = data[2], times = data[3], lambda1 = data[4], lambda2 = data[5], thrh = data[6])
numT = len(data.times)

chromLen = 3000000
chromNum = 1000
binSize = 100

theta = 2*chromLen*data.thrh[0]/binSize
rho = 2*chromLen*data.thrh[1]/binSize

mscl = " 4 " + str(chromNum) + " -t " + str(theta) + " -r " + str(rho) + " " + str(chromLen) + " -l -I 2 2 2 "
for i in range(data.splitT):
    mscl += " -en " + str(data.times[i]/2.0) + " 1 " + str(1.0/data.lambda1[i])
    mscl += " -en " + str(data.times[i]/2.0) + " 2 " + str(1.0/data.lambda2[i])

mscl += " -em " + str(data.times[ data.migStart ]/2.0) + " 1 2 " + str(2*data.mu[0])
mscl += " -em " + str(data.times[ data.migStart ]/2.0) + " 2 1 " + str(2*data.mu[1])
mscl += " -eM " + str(data.times[ data.migEnd ]/2.0) + " 0.0 "

mscl += " -ej " + str(data.times[ data.splitT ]/2.0) + " 2 1 "
mscl += " -eM " + str(data.times[ data.splitT ]/2.0) + " 0.0 "

for i in range(data.splitT, numT):
    mscl += " -eN " + str(data.times[i]/2.0) + " " + str(1.0/data.lambda1[i])

print(mscl)

#migrationIO.PlotMS(fms)

#./foreign/msHOT-lite/msHOT-lite 4 1000 -t 8196 -r 1355 3000000 -l -I 2 2 2 -n 2 1.0 -em 0.0 1 2 2.0 -em 0.0 2 1 2.0 -en 0.01 1 0.05 -en 0.01 2 0.05 -en 0.0375 1 0.5 -en 0.0375 2 0.5 -ej 1.25 2 1 -eM 1.25 0.0 -eN 1.25 1.0 