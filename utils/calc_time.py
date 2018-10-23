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
cur_dir = os.path.dirname(os.path.abspath(__file__))
mig_dir = os.path.join(cur_dir, '..')
sys.path.append(mig_dir)
import os
import collections
import argparse
import numpy
import math
from math import (exp,log)
import migrationIO

parser = argparse.ArgumentParser(description='Migration inference from PSMC.')

parser.add_argument('fpsmc1',
                    help='psmc file 1')
parser.add_argument('fpsmc2',
                    help='psmc file 2')

parser.add_argument('-wd', nargs=1, default='',
                    help='working directory (path to data files)')
                    
parser.add_argument('--sdate', nargs=1, type=float, default=0,
                    help='dating of the second sample (for ancient genome, in years - units should be set properly)') 
parser.add_argument('-rd', nargs=1, type=int, default=-1,
                    help='Round (RD) in PSMC file (default -1 for the last round, in this case the number of rounds should be exactly the same in both files)')
parser.add_argument('--funits', nargs=1, type=str, default="setunits.txt",
                    help='File name with units to be used to rescale times and EPS.')

clargs = parser.parse_args()
if isinstance(clargs.wd, list):
    clargs.wd = clargs.wd[0]
if isinstance(clargs.rd, list):
    clargs.rd = clargs.rd[0]
if isinstance(clargs.sdate, list):
    clargs.sdate = clargs.sdate[0]
if isinstance(clargs.funits, list):
    clargs.funits = clargs.funits[0]
    
units = migrationIO.Units()
units.SetUnitsFromFile(clargs.funits)

fpsmc1 = os.path.join( clargs.wd, clargs.fpsmc1 )
fpsmc2 = os.path.join( clargs.wd, clargs.fpsmc2 )

inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2, clargs.sdate, clargs.rd, False)
for splitT in range(len(inputData[0])):
    print("splitT = ", splitT, "\ttime = ", int(sum(inputData[0][0:int(splitT)])*inputData[2]))
