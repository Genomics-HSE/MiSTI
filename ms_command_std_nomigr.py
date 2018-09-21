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
hapLen = 3*1000*1000
chromNum = 1000
#chromNum = 300

path = "./foreign/msHOT-lite/"
simul = "msHOT-lite"
simul = path+simul+" "
simul = ""
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
N = [N_0, 0.05*N_0, 2.0*N_0, 0.5*N_0]
#Migration
mu = [0.00005, 0.000075]
#mu = [0.0, 0.0]
#Split times
T = [4.0*N_0*0.01, 4.0*N_0*0.0375, 4.0*N_0*1.25]


th = 6.83*10**(-8) #Mutation rate per site per generation
rh = 1.13*10**(-8) #Recombination rate per site per generation
#6.83e-8

for i in range(len(N)):
    N[i] = N[i]/N_0
    
for i in range(len(T)):
    T[i] = T[i]/(4.0*N_0)

print 'MS command generator for two populations with migration'
    
command = '' + simul + `int(popNum)` + ' ' + `int(chromNum)` + ' '

theta = int(4*th*N_0*hapLen)
rho = int(4*rh*N_0*hapLen)

command += '-t ' + `theta` + ' -r ' + `rho` + ' ' + `hapLen` + ' ' + additionalParams + ' '

command += '-I 2 ' + `pop[0]` + ' ' + `pop[1]` + ' '
command += '-n 2 ' + `N[0]` + ' '
mu[0] = 0.0#0.5, 1.0, 2.0
mu[1] = 0.0
command += '-em 0.0 1 2 ' + `mu[0]` + ' '# + `4*N_0*mu[0]` + ' '
command += '-em 0.0 2 1 ' + `mu[1]` + ' '# + `4*N_0*mu[1]` + ' '

#command += '-eM ' + `T[0]` + ' 0.0 '
command += '-en ' + `T[0]` + ' 1 ' + `N[1]` + ' '
command += '-en ' + `T[0]` + ' 2 ' + `N[1]` + ' '

command += '-en ' + `T[1]` + ' 1 ' + `N[3]` + ' '
command += '-en ' + `T[1]` + ' 2 ' + `N[3]` + ' '

command += '-ej ' + `T[2]` + ' 2 1 '
command += '-eM ' + `T[2]` + ' 0.0 '
command += '-eN ' + `T[2]` + ' ' + `N[0]` + ' '
#command += '> data_std_migr_same/sim.ms'
m0 = `mu[0]`
m1 = `mu[1]`
m0 = m0.replace(".","")
m1 = m1.replace(".","")
path = "data_" + "bn" + m0 + "_" + "bn" + m1
command = "./run_sim.sh " + path + " " + "\"" + command + "\""
print command

sys.exit(1)
print "\n\nMigration inference parameters"
lap = 'lambdas = ['
for i in range(3):
    lap += '[' + `1.0/N[2*i]` + ', ' + `1.0/N[2*i+1]` + ']'
    if i != 2:
        lap += ', '
lap += ']'
print lap
print 'mu = [' + `mu[0]*N_0` + ', ' + `mu[1]*N_0` + ']'
print 'times = [' + `T[0]*2` + ', '  + `T[1]*2` + ']'m