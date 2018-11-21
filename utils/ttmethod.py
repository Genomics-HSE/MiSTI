#!/usr/bin/env python3

#tt-method is described in the paper "Southern African ancient genomes estimate modern human divergence to 350,000 to 260,000 years ago" by Schlebusch et al, Genetics 2017
#implemented by Vladimir Shchur vlshchur@gmail.com
#input file format is consistent with joint allele frequency spectrum format used in MiSTI (https://github.com/vlshchur/MiSTI)
#haplotype size is specified by user

import sys
import os
import argparse
import numpy
from migrationIO import ReadJAFS
from math import log

parser = argparse.ArgumentParser(description='Implementation of TT-method (Schlebusch et al, Genetics 2017)')

parser.add_argument('jafs',
                    help='Joint allele frequency spectrum')
parser.add_argument('haplen',
                    help='Halotype length (total number of sites, both variable and non-variable)')

parser.add_argument('-y', nargs=1, type=float, default=1,
                    help='years per generation')
parser.add_argument('-mu', nargs=1, type=float, default=1.25e-8,
                    help='mutation rate (per basepair per generation)')

clargs = parser.parse_args()

if isinstance(clargs.y, list):
    clargs.y = clargs.y[0]
if isinstance(clargs.mu, list):
    clargs.mu = clargs.mu[0]

spectrum = ReadJAFS(clargs.jafs).jafs
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

a1 = 2*m5/(2*m6+m5)
a2 = 2*m5/(2*m7+m5)

theta = 3/M*(2*m6+m5)*(2*m7+m5)/(8*m5)/2
theta1 = -T1/log(a1)/2
theta2 = -T2/log(a2)/2

mu = clargs.mu
Ygen = clargs.y
T1y = T1/mu*Ygen
T2y = T2/mu*Ygen

N_A = theta/mu
N_1 = theta1/mu
N_2 = theta2/mu

print("Implementation of tt method (Schlebusch et al, Genetics 2017)")
print("T1 = ", T1y)
print("T2 = ", T2y)
print("N_A = ", N_A, "\tN_1 = ", N_1, "\tN_2 = ", N_2)

