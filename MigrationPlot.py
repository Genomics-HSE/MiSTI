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

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('fpsmc1',
                    help='psmc file 1')
parser.add_argument('fpsmc2',
                    help='psmc file 2')
parser.add_argument('fmigr',
                    help='migr file')

parser.add_argument('-wd', nargs=1, default='',
                    help='working directory (path to data files)')
parser.add_argument('-o', nargs=1, default='plot.pdf',
                    help='output filename')

parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')

clargs = parser.parse_args()
if isinstance(clargs.wd, list):
    clargs.wd = clargs.wd[0]
if isinstance(clargs.o, list):
    clargs.o = clargs.o[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]

RD = clargs.rd
fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc1 )
fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc2 )
fmigr  = os.path.join( clargs.wd, clargs.fmigr  )
fout   = os.path.join( clargs.wd, clargs.o )

print("Output file: ", fout)

migrationIO.PlotInit()

data = migrationIO.ReadPSMC(fpsmc1, fpsmc2, RD, True)
migrationIO.ReadMigration(fmigr, True, data[2], data[3])
#migrationIO.PlotMS(fms)

migrationIO.SavePlot(fout)