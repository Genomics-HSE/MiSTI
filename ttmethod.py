#!/Users/shchur/anaconda3/bin/python3

#!/usr/bin/env python3

import sys
import os
import collections
import argparse
import numpy
import math
from migrationIO import ReadJAFS


parser = argparse.ArgumentParser(description='Implementation of TT-method (Schlebusch et al, Genetics 2017)')

parser.add_argument('jafs',
                    help='Joint allele frequency spectrum')
parser.add_argument('haplen',
                    help='Halotype length (total number of sites, both variable and non-variable)')

clargs = parser.parse_args()

spectrum = ReadJAFS(clargs.jafs)

#Converting to paper notations
M = clargs.haplen
m1 = spectrum[0]
m2 = spectrum[2]
m3 = spectrum[1]
m4 = spectrum[5]
m5 = spectrum[3]
m6 = spectrum[4]
m7 = spectrum[6]

T1 = (m1/2+m3-(2*m6+m5)*(6*m7+m5)/8/m5)/M
T2 = (m2/2+m4-(2*m7+m5)*(6*m6+m5)/8/m5)/M

print("T1 = ", T1)
print("T2 = ", T2)