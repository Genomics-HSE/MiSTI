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

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('fmigr',
                    help='migr file')

parser.add_argument('--fpsmc', '-fp', nargs=2, type=str, default=None,
                    help='psmc files')

parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')

parser.add_argument('-wd', nargs=1, default='',
                    help='working directory (path to data files)')
parser.add_argument('-o', nargs=1, default='plot.pdf',
                    help='output filename')
                    
parser.add_argument('--sdate', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome, in years - units should be set properly)') 
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('--maxY', nargs=1, type=float, default=None,
                    help='Range for Y axis (upper bound).')
parser.add_argument('--minY', nargs=1, type=float, default=None,
                    help='Range for Y axis (upper bound).')
parser.add_argument('--maxX', nargs=1, type=float, default=None,
                    help='Range for Y axis (upper bound).')
parser.add_argument('--minX', nargs=1, type=float, default=None,
                    help='Range for Y axis (upper bound).')

parser.add_argument('--hideProbs', '-hp', action='store_true',
                    help='Hide probability pannels from the plot.')
                    
clargs = parser.parse_args()
if isinstance(clargs.wd, list):
    clargs.wd = clargs.wd[0]
if isinstance(clargs.o, list):
    clargs.o = clargs.o[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]
if isinstance(clargs.maxY, list):
    clargs.maxY = clargs.maxY[0]
if isinstance(clargs.minY, list):
    clargs.minY = clargs.minY[0]
if isinstance(clargs.maxX, list):
    clargs.maxX = clargs.maxX[0]
if isinstance(clargs.minX, list):
    clargs.minX = clargs.minX[0]
if isinstance(clargs.maxY, list):
    clargs.maxY = clargs.maxY[0]
if isinstance(clargs.sdate, list):
    clargs.sdate = clargs.sdate[0]
if isinstance(clargs.hideProbs, list):
    clargs.hideProbs = clargs.hideProbs[0]
    
if isinstance(clargs.trueEPS, list):
    clargs.trueEPS = clargs.trueEPS[0]
    
units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)
units.PrintUnits()

if clargs.fpsmc is not None:
    fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc[0] )
    fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc[1] )
fmigr  = os.path.join( clargs.wd, clargs.fmigr  )
fout   = os.path.join( clargs.wd, clargs.o )

print("Output file: ", fout)

migrationIO.PlotInit(hideProbs = clargs.hideProbs)

if clargs.fpsmc is not None:
    data = migrationIO.ReadPSMC(fpsmc1, fpsmc2, clargs.sdate, clargs.rd, True, maxY = clargs.maxY)
else:
    data = [None, None, 1, 1]
plotLimits = {
    "maxY": clargs.maxY,
    "minY": largs.minY,
    "maxX": clargs.maxX,
    "minX": clargs.minX
}
migrationIO.ReadMigration(fmigr, True, data[2], data[3])
#migrationIO.PlotMS(fms)
#plt.legend()
migrationIO.SavePlot(fout, plotLimits)
