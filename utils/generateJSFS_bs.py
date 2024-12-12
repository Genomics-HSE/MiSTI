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
from migrationIO import PrintJAFSFile, ReadJAFS, BootstrapJAFS

if len(sys.argv) != 3:
    print("./generateJSFS_bs.py <number of bs samples> <Joint SFS file with multiple chunks>")
    exit(0)

bs_size = int(sys.argv[1])
fjafs = sys.argv[2]

pop1, pop2 = [], []

dataJAFS = ReadJAFS(sys.argv[2], True)

trueSFS = [0 for _ in range(8)]
for sfs in dataJAFS.jafs:
    trueSFS = [v+u for v, u in zip(trueSFS, sfs)]

jafs = [trueSFS]

for i in range(bs_size):
    jafs.append( BootstrapJAFS(dataJAFS) )

PrintJAFSFile(jafs, dataJAFS.pop1, dataJAFS.pop2)
