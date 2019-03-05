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

pop1, pop2 = [], []


jafs = []
for fn in sys.argv[1:]:
    if os.path.isdir(fn):
        for fn1 in os.listdir(fn):
            if fn1[0] != ".":
                dataJAFS = ReadJAFS(os.path.join( fn, fn1 ), True)
                jaf = dataJAFS.jafs
                jafs.append(jaf)
                if dataJAFS.pop1 is not None:
                    pop1.append(dataJAFS.pop1)
                if dataJAFS.pop2 is not None:
                    pop2.append(dataJAFS.pop2)
    else:
        dataJAFS = ReadJAFS(fn, True)
        jaf = dataJAFS.jafs
        jafs.append(jaf)
        if dataJAFS.pop1 is not None:
            pop1.append(dataJAFS.pop1)
        if dataJAFS.pop2 is not None:
            pop2.append(dataJAFS.pop2)
pop1 = list(set(pop1))
pop2 = list(set(pop2))
pop1s = "+".join(pop1)
pop2s = "+".join(pop2)

PrintJAFSFile(jaf, pop1s, pop2s)
