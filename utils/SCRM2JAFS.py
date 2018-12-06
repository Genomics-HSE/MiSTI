#!/usr/bin/env python3

#    Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)
#This script produces many jafs from one-to-one comparison of individuals.

import sys
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
mig_dir = os.path.join(cur_dir, '..')
sys.path.append(mig_dir)
from migrationIO import PrintJAFSFile


if len(sys.argv) < 2:
    print("./SCRM2JAFS.py <INPUT FILE>")
    exit(0)
hapNum = 40
chrLim = 1000


fn = sys.argv[1]
sfs = [0 for _ in range(hapNum)]
chrs = [[None] for _ in range(hapNum)]
chrPr = 0
jaf = [0 for _ in range(7)]

with open(fn) as f:
    for line in f:
        if line[0:9] == "positions":
            for i in range(4):
                line = next(f)
                chrs[i] = [int(v) for v in line[0:-1]]
    
            for i in range(len(chrs[0])):
                chrs[0][i]
                s0 = int( chrs[0][i] ) + int( chrs[1][i] )
                s1 = int( chrs[2][i] ) + int( chrs[3][i] )
    
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

total = sum(jaf)
for v in jaf:
    print(v/total)
#PrintJAFSFile(jaf)
