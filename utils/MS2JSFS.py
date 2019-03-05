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




def AddPosition(jaf, jafs, chLen, begin, end, remaining_ch_len):
    if end - begin < remaining_ch_len:
        remaining_ch_len -= (end - begin)
    else:
        SaveJAF(jaf, jafs, chLen)
#        chunks_processed += 1
        tmp_len =  (end - begin) - remaining_ch_len
        remaining_ch_len = chunkLen
        remaining_ch_len -= tmp_len
    return remaining_ch_len

def SaveJAF(jaf, jafs, chLen):
    jafs.append([v for v in jaf])
    jafs[-1].insert(0, chLen)
    for i in range(len(jaf)):
        jaf[i] = 0

parser = argparse.ArgumentParser(description='This script calculates joint SFS from Heng Li\'s output format of msHOT-lite (-l option).')

parser.add_argument('inputfile',
                    help='input file (msHOT-lite with -l option)')
parser.add_argument('-p', nargs=2, type=str,
                    help='population names (makes it easier to check the consistence of the pipeline.)')
parser.add_argument('-n', nargs=1, type=int, default=200,
                    help='Number of chunks for bootstrap.')


clargs = parser.parse_args()

if isinstance(clargs.n, list):
    clargs.n = clargs.n[0]

#if len(sys.argv) < 2:
#    print("./MS2JAF.py <INPUT FILE> [pop1] [pop2]")
#    exit(0)

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
jafs = []

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
    chunkLen = math.ceil(numChrom*chromLen/clargs.n)
    pr_position = 0
    ch_len = chunkLen
#    chunks_processed = 0
    while line != "EOF":
        while not (line.startswith("@begin") or line == "EOF"):
            line = next(f, "EOF")
            if line.startswith("segsites:"):
                ch_len = AddPosition(jaf, jafs, chunkLen, 0, chromLen, ch_len)
#                if chromLen < ch_len:
#                    ch_len -= chromLen
#                else:
#                    SaveJAF(jaf, jafs, chromLen)
#                    chunks_processed += 1
#                    tmp_len =  chromLen - ch_len
#                    ch_len = chunkLen
#                    ch_len -= tmp_len
        while not (line.startswith("@end") or line == "EOF"):
            line = next(f, "EOF").rstrip("\n")
            pars = line.split("\t")
            if line.startswith("@end"):
                ch_len = AddPosition(jaf, jafs, chunkLen, pr_position, chromLen, ch_len)
#                if position - pr_position < ch_len:
#                    ch_len -= (position - pr_position)
#                else:
#                    SaveJAF(jaf, jafs, chromLen)
#                    chunks_processed += 1
#                    tmp_len =  (position - pr_position) - ch_len
#                    ch_len = chunkLen
#                    ch_len -= tmp_len
                pr_position = 0
            if len(pars) != 2:
                continue
            position = int(pars[0])
            ch_len = AddPosition(jaf, jafs, chunkLen, pr_position, position, ch_len)
#            if position - pr_position < ch_len:
#                ch_len -= (position - pr_position)
#            else:
#                SaveJAF(jaf, jafs, chromLen)
#                chunks_processed += 1
#                tmp_len =  (position - pr_position) - ch_len
#                ch_len = chunkLen
#                ch_len -= tmp_len
#                pr_position = position
            pr_position = position
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
#print("chunks_processed", chunks_processed)
#if chunks_processed == clargs.n:
#    print("IS IT EXPECTED?? CHECK ME!!!")
#    sys.exit(0)
if len(jafs) != clargs.n:
    SaveJAF(jaf, jafs, chunkLen - ch_len)
PrintJAFSFile(jafs, pop1, pop2)