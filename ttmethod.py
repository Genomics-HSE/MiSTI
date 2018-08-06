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

parser.add_argument('-y', nargs=1, type=float, default=30,
                    help='years per generation')
parser.add_argument('-mu', nargs=1, type=float, default=1.25e-8,
                    help='mutation rate (per basepair per generation)')

clargs = parser.parse_args()

if isinstance(clargs.y, list):
    clargs.y = clargs.y[0]
if isinstance(clargs.fout, list):
    clargs.mu = clargs.mu[0]

spectrum = ReadJAFS(clargs.jafs)
spectrum = spectrum[1:]
#Converting to paper notations
M = float(clargs.haplen)
m1 = spectrum[0]
m2 = spectrum[2]
m3 = spectrum[1]
m4 = spectrum[5]
m5 = spectrum[3]
m6 = spectrum[4]
m7 = spectrum[6]

T1 = (m1/2+m3-(2*m6+m5)*(6*m7+m5)/8/m5)/M
T2 = (m2/2+m4-(2*m7+m5)*(6*m6+m5)/8/m5)/M

mu = clargs.mu
Ygen = clargs.y
T1y = T1/mu*Ygen
T2y = T2/mu*Ygen

print("T1 = ", T1y)
print("T2 = ", T2y)
