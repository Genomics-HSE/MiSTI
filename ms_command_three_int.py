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
hapLen = 30*1000*100
chromNum = 1000

simul = "msHOT-lite"
additionalParams = "-l"

for i in range(1, len(sys.argv)):
    if sys.argv[i] == "-scrm":
        simul = "scrm"
        linkage = "-l 100"
    if sys.argv[i] == "-ms":
        additionalParams = "ms"
        additionalParams = ""
    if sys.argv[i] == "-chr":
        i += 1
        chromNum = int(sys.argv[i])

#Effective population sizes
N_0 = 10000.0
N = [N_0, 4000.0, 7000.0, 6000.0, 5000.0, 5000.0]
#Migration
mu = [0.01, 0.003]
#mu = [0.0, 0.0]
#Split times
T = [5000.0, 12000.0]


th = 10**(-7) #Mutation rate per site per generation
rh = 10**(-8) #Recombination rate per site per generation
th = 6.83*10**(-8)
rh = 1.13*10**(-8)

for i in range(len(N)):
    N[i] = N[i]/N_0
    
for i in range(len(T)):
    T[i] = T[i]/(4.0*N_0)

print 'MS command generator for two populations with migration'
    
command = '' + simul + ' ' + `int(popNum)` + ' ' + `int(chromNum)` + ' '

theta = int(4*th*N_0*hapLen)
rho = int(4*rh*N_0*hapLen)

command += '-t ' + `theta` + ' -r ' + `rho` + ' ' + `hapLen` + ' ' + additionalParams + ' '

command += '-I 2 ' + `pop[0]` + ' ' + `pop[1]` + ' '
command += '-n 2 ' + `N[1]` + ' '
command += '-em 0.0 1 2 ' + `4*N_0*mu[0]` + ' '
command += '-em 0.0 2 1 ' + `4*N_0*mu[1]` + ' '

#command += '-eM ' + `T[0]` + ' 0.0 '
command += '-en ' + `T[0]` + ' 1 ' + `N[2]` + ' '
command += '-en ' + `T[0]` + ' 2 ' + `N[3]` + ' '

command += '-ej ' + `T[1]` + ' 2 1 '
command += '-eN ' + `T[1]` + ' ' + `N[4]` + ' '
command += '> data/sim_3_migr.ms'
print command


print "\n\nMigration inference parameters"
lap = 'lambdas = ['
for i in range(3):
    lap += '[' + `1.0/N[2*i]` + ', ' + `1.0/N[2*i+1]` + ']'
    if i != 2:
        lap += ', '
lap += ']'
print lap
print 'mu = [' + `mu[0]*N_0` + ', ' + `mu[1]*N_0` + ']'
print 'times = [' + `T[0]*2` + ', '  + `T[1]*2` + ']'