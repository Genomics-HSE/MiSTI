#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    print("./MS2JAF.py <INPUT FILE>")
    exit(0)

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
            fr = [int(v) for v in line[1]]
            
            s0 = int(fr[0]) + int(fr[1])
            s1 = int(fr[2]) + int(fr[3])
            
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

norm = sum(jaf)
jaf = [v/norm for v in jaf]
print(jaf)