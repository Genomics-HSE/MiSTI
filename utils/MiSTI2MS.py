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
scaling = migrationIO.SetScaling()
#data = MigData(splitT = data[0], migStart = data[1], migEnd = data[2], times = data[3], lambda1 = data[4], lambda2 = data[5], thrh = data[6])
numT = len(data.times)

chromLen = 3000000
chromNum = 1000
binSize = 100


N0 = data.thrh[0]/(4*scaling[0]*scaling[1])
N0_rescale = 10000/N0
theta = chromLen*data.thrh[0]/binSize*N0_rescale
rho = chromLen*data.thrh[1]/binSize*N0_rescale

mscl = " 4 " + str(chromNum) + " -t " + str(theta) + " -r " + str(rho) + " " + str(chromLen) + " -l -I 2 2 2 "
lp = [0, 0]
for i in range(data.splitT):
    if lp[0] != data.lambda1[i]:
        mscl += " -en " + str(data.times[i]/2.0/N0_rescale) + " 1 " + str(N0_rescale/data.lambda1[i])
        lp[0] = data.lambda1[i]
    if lp[1] != data.lambda2[i]:
        mscl += " -en " + str(data.times[i]/2.0/N0_rescale) + " 2 " + str(N0_rescale/data.lambda2[i])
        lp[1] = data.lambda2[i]

mscl += " -em " + str(data.times[ data.migStart ]/2.0/N0_rescale) + " 1 2 " + str(2*data.mu[0]*N0_rescale)
mscl += " -em " + str(data.times[ data.migStart ]/2.0/N0_rescale) + " 2 1 " + str(2*data.mu[1]*N0_rescale)
mscl += " -eM " + str(data.times[ data.migEnd ]/2.0/N0_rescale) + " 0.0 "

mscl += " -ej " + str(data.times[ data.splitT ]/2.0/N0_rescale) + " 2 1 "
mscl += " -eM " + str(data.times[ data.splitT ]/2.0/N0_rescale) + " 0.0 "

lp = 0
for i in range(data.splitT, numT):
    if lp != data.lambda1[i]:
        mscl += " -eN " + str(data.times[i]/2.0/N0_rescale) + " " + str(N0_rescale/data.lambda1[i])
        lp = data.lambda1[i]

print(mscl)

#migrationIO.PlotMS(fms)

#./foreign/msHOT-lite/msHOT-lite 4 1000 -t 8196 -r 1355 3000000 -l -I 2 2 2 -n 2 1.0 -em 0.0 1 2 2.0 -em 0.0 2 1 2.0 -en 0.01 1 0.05 -en 0.01 2 0.05 -en 0.0375 1 0.5 -en 0.0375 2 0.5 -ej 1.25 2 1 -eM 1.25 0.0 -eN 1.25 1.0 