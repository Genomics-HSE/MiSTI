#!/usr/bin/env python3
import sys
import os
cur_dir = os.path.dirname(os.path.abspath(__file__))
mig_dir = os.path.join(cur_dir, '..')
sys.path.append(mig_dir)
from migrationIO import PrintJAFSFile, PrintErr

if len(sys.argv) < 2:
    print("./ANGSDSFS.py <INPUT FILE> [pop1 pop2]")
    exit(0)
    
pop1, pop2 = False, False

if len(sys.argv) == 4:
    pop1 = sys.argv[2]
    pop2 = sys.argv[3]
else:
    PrintErr("IMPORTANT NOTICE!!! It is strongly recommended to supply population 1 and population 2 names to ensure that the order of psmc files is not swapped relatively to the joint allele frequency spectrum.")

fn = sys.argv[1]
jafs = []

with open(fn) as f:
    for line in f:
        line = line.rstrip("\n")
        freqs = line.split(" ")
        sfs = [float(v) for v in freqs[0:8]]
        jafs.append([
            sum(sfs),
            sfs[3],
            sfs[6],
            sfs[1],
            sfs[4],
            sfs[7],
            sfs[2],
            sfs[5]
        ])
PrintJAFSFile(jafs, pop1, pop2)
