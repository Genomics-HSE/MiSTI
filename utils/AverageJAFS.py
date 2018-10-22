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
from migrationIO import PrintJAFSFile, ReadJAFS

if len(sys.argv) < 2:
    print("./MS2JAF.py <ANY NUMBER OF INPUT FILES OR DIRECTORIES>")
    exit(0)

pop1, pop2 = False, False

jaf = [0 for _ in range(7)]
tf = 0
for fn in sys.argv[1:]:
    if os.path.isdir(fn):
        for fn1 in os.listdir(fn):
            if fn1[0] != ".":
                jafs = ReadJAFS( os.path.join( fn, fn1 ) )[1:]
                jaf = [u + v for u, v in zip(jaf, jafs)]
                tf += 1
    else:
        jafs = ReadJAFS(fn)[1:]
        jaf = [u + v for u, v in zip(jaf, jafs)]
        tf += 1
jaf = [int(v/tf) for v in jaf]

PrintJAFSFile(jaf, pop1, pop2)
