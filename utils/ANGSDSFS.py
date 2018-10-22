#!/usr/bin/env python3
import sys
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
jaf = [0 for _ in range(7)]

with open(fn) as f:
    for line in f:
        freqs = line.split(" ")
        for i in range(1, 8):
            jaf[i-1] += float( freqs[i] )

jaf = [round(u) for u in jaf]
jaf = [
    jaf[2],
    jaf[5],
    jaf[0],
    jaf[3],
    jaf[6],
    jaf[1],
    jaf[4]
]
PrintJAFSFile(jaf, pop1, pop2)
