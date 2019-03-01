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
import argparse
import math
cur_dir = os.path.dirname(os.path.abspath(__file__))
mig_dir = os.path.join(cur_dir, '..')
sys.path.append(mig_dir)
from migrationIO import PrintJAFSFile, PrintErr



parser = argparse.ArgumentParser(description='This script calculates joint SFS from Heng Li\'s output format of msHOT-lite (-l option).')

parser.add_argument('inputfile',
                    help='input file (msHOT-lite with -l option)')
parser.add_argument('-p', nargs=2, type=str,
                    help='population names (makes it easier to check the consistence of the pipeline.)')
parser.add_argument('-n', nargs=1, type=int, default=200,
                    help='File name with units to be used to rescale times and EPS.')


clargs = parser.parse_args()

if len(sys.argv) < 2:
    print("./MS2JAF.py <INPUT FILE> [pop1] [pop2]")
    exit(0)

pop1, pop2 = False, False

if clargs.p is not None:
    pop1 = clargs.p[0]
    pop2 = clargs.p[1]

#Can be used for a permutation of haplotypes
h0 = 0
h1 = 1
h2 = 2
h3 = 3

fn = sys.argv[1]
jaf = [0 for _ in range(7)]

with open(fn) as f:
    line = next(f, "EOF")
    if line != "EOF":
        pars = line.split(" ")
        numChrom = int(pars[2])
        chromLen = 0
        for i in range(len(pars)):
            if pars[i] == "-r" and i+2 < len(pars):
                chromLen = int(pars[i+2])
                break
        if chromLen <= 0:
            PrintErr("Unknown number of chromosomes. The script is designed to work with ms commands containing -r argument.")
            sys.exit(0)
    chunkLen = math.floor(numChrom*chromLen/clargs.n)
    while line != "EOF":
        while not (line.startswith("@begin") or line == "EOF"):
            line = next(f, "EOF")
            if line.startswith("segsites:"):
                1+1
        while not (line.startswith("@end") or line == "EOF"):
            line = next(f, "EOF").rstrip("\n")
            pars = line.split("\t")
            if len(pars) != 2:
                continue
            fr = list(pars[1][0:4])
            
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