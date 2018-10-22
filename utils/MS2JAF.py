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
from migrationIO import PrintJAFSFile

if len(sys.argv) < 2:
    print("./MS2JAF.py <INPUT FILE> [pop1] [pop2]")
    exit(0)

pop1, pop2 = False, False

if len(sys.argv) == 4:
    pop1 = sys.argv[2]
    pop2 = sys.argv[3]

h0 = 0
h1 = 1
h2 = 2
h3 = 3
if len(sys.argv) == 6:
    h0 = int(sys.argv[2])
    h1 = int(sys.argv[3])
    h2 = int(sys.argv[4])
    h3 = int(sys.argv[5])
    if h0 + h1 + h2 + h3 != 6 or min(h0, h1, h2, h3) < 0:
        print("Haplotype numbers should be a permutation of 0, 1, 2 and 3.")
        sys.exit(0)

fn = sys.argv[1]
jaf = [0 for _ in range(7)]

with open(fn) as f:
    for line in f:        
        for _ in range(2):
            line = next(f)
        chrLen = int( next(f) )
        lc = 0
        while 1:
            line = next(f)
            lc += 1
            if lc > int(chrLen):
                print("Too many segsites, expected at most " + str(chrLen))
                sys.exit(0)
            if line == "@end\n":
                break
            line = line.split("\t")
            fr = list(line[1][0:4])
            
            s0 = int(fr[h0]) + int(fr[h1])
            s1 = int(fr[h2]) + int(fr[h3])
            
            if s0 == 0:
                if s1 == 1:
                    jaf[2] += 1
                elif s1 == 2:
                    jaf[5] += 1
            if s0 == 1:
                if s1 == 0:
                    jaf[0] += 1
                elif s1 == 1:
                    jaf[3] += 1
                elif s1 == 2:
                    jaf[6] += 1
            if s0 == 2:
                if s1 == 0:
                    jaf[1] += 1
                elif s1 == 1:
                    jaf[4] += 1

PrintJAFSFile(jaf, pop1, pop2)