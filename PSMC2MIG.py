#!/usr/bin/env python3
import sys

if len(sys.argv) < 2:
    print("./PSMC2MIG.py <INPUT FILENAME> <RD - optional>")
    print("See PSMC format for explanation of RD.")
    sys.exit(0)

fn = sys.argv[1]
RD = -1
maxRD = -1
Tk = []
Lk = []
if len(sys.argv) == 3:
     RD = int( sys.argv[2] )

with open(fn) as f:
    for line in f:
        line = line.split()
        if line[0] == "RD":
            maxRD = int( line[1] )
    
    if maxRD == -1:
        print("Corrupted of empty input file")
        sys.exit(0)

    if RD == -1:
        RD = maxRD

with open(fn) as f:  
    for line in f:
        line = line.split()
        if line[0] != "RD" or int( line[1] ) != RD:
            continue
        while line[0] != "RS":
            line = next(f)
            line = line.split()
        while line[0] != "PA":
            if line[0] != "RS":
                print("Unexpected line.")
                sys.exit(0)
            Tk.append( float(line[2]) )
            Lk.append( float(line[3]) )
            line = next(f)
            line = line.split()
        break

print(Tk)
print(Lk)
