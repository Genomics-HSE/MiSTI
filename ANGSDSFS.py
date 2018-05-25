#!/usr/bin/env python3
import sys
import migrationIO

if len(sys.argv) < 2:
    print("./ANGSDSFS.py <INPUT FILE>")
    exit(0)

fn = sys.argv[1]
jaf = [0 for _ in range(7)]

with open(fn) as f:
    for line in f:
        freqs = line.split("\t")
        for i in range(1, 9):
            jaf[i-1] += float( freqs[i] )

PrintJAFSFile(jaf)
