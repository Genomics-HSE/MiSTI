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

cur_dir = os.path.dirname(os.path.abspath(__file__))
mig_dir = os.path.join(cur_dir, '..')
sys.path.append(mig_dir)
import migrationIO
plt.ioff()

parser = argparse.ArgumentParser(description='Migration inference from PSMC help scripts.')
parser.add_argument('ms',
                    help='ms command line')
clargs = parser.parse_args()

inputData = migrationIO.ReadMS(clargs.ms)

print(inputData)
