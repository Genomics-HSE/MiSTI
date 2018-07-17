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
import MigData from migrationIO
plt.ioff()

parser = argparse.ArgumentParser(description='Migration inference from PSMC help scripts.')
parser.add_argument('fmigr',
                    help='migr file')
clargs = parser.parse_args()
fmigr  = os.path.join( clargs.wd, clargs.fmigr  )

data = migrationIO.ReadMigration(fmigr, False)
data = MigData(splitT = data[0], migStart = data[1], migEnd = data[2], times = data[3], lambda1 = data[4], lambda2 = data[5], thrh = data[6])
numT = len(data.times)

chromLen = 30000000
chromNum = 100
binSize = 100

theta = chromLen*data.thrh[0]/binSize

theta = 
mscl = ""
for i in range(numT):
    mscl += " -en " + str(data.times[i]) + " 1 " + str(data.lambda1[i])
    mscl += " -en " + str(data.times[i]) + " 2 " + str(data.lambda2[i])

#migrationIO.PlotMS(fms)




command += '-em 0.0 1 2 ' + `mu[0]` + ' '# + `4*N_0*mu[0]` + ' '
command += '-em 0.0 2 1 ' + `mu[1]` + ' '# + `4*N_0*mu[1]` + ' '

#command += '-eM ' + `T[0]` + ' 0.0 '
command += '-en ' + `T[0]` + ' 1 ' + `N[1]` + ' '
command += '-en ' + `T[0]` + ' 2 ' + `N[1]` + ' '

command += '-en ' + `T[1]` + ' 1 ' + `N[3]` + ' '
command += '-en ' + `T[1]` + ' 2 ' + `N[3]` + ' '

command += '-ej ' + `T[2]` + ' 2 1 '
command += '-eM ' + `T[2]` + ' 0.0 '