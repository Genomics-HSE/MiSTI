#!/usr/bin/python
import sys
'''
import re
import os
import csv
import math
import numpy
import random'''

'''popNum = False
hapLen = False
N0 = False
N1 = False
pop1Num = False'''

#Population sample sizes
pop = [2, 2]
popNum = sum(pop)
hapLen = 1000*2500
chromNum = 5

simul = "ms"
linkage = ""

for i in range(1, len(sys.argv)):
    if sys.argv[i] == "-scrm":
        simul = "scrm"
        linkage = "-l 100 "
    if sys.argv[i] == "-chr":
        i += 1
        chromNum = int(sys.argv[i])

#Effective population sizes
N_0 = 10000.0
N = [N_0, 4000.0]
#Migration
mu = [0.01, 0.003]
#Split times
T = [5000.0]


th = 10**(-8) #Mutation rate per site per generation
rh = 10**(-8) #Recombination rate per site per generation

for i in range(len(N)):
    N[i] = N[i]/N_0
    
for i in range(len(T)):
    T[i] = T[i]/(4.0*N_0)

print 'MS command generator for two populations with migration'
    
command = './' + simul + ' ' + `int(popNum)` + ' ' + `int(chromNum)` + ' '

theta = 4*th*N_0*hapLen
rho = 4*rh*N_0*hapLen

command += '-t ' + `theta` + ' -r ' + `rho` + ' ' + `hapLen` + ' ' + linkage

command += '-eN ' + `T[0]` + ' ' + `N[1]` + ' '
command += '> sim.txt'
print command
print 'lambdas = [[' + `1.0/N[0]` + ', ' + `1.0/N[0]` + ']' + ", " + '[' + `1.0/N[1]` + ', ' + `1.0/N[1]` + ']]'
print 'mu = [' + `mu[0]` + ', ' + `mu[1]` + ']'
print 'times = [' + `T[0]*2` + ']'