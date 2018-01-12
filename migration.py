#!/usr/bin/env python3

import sys
import collections
import numpy
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)

from CorrectLambda import CorrectLambda
from TwoPopulations import TwoPopulations
from OnePopulation import OnePopulation
import migrationIO

class MigrationInference:
    def __init__(self, times, lambdas, dataJAFS, mu, splitT, theta, correct = True):
        self.correct = correct
        #Model parameters
        self.theta = theta#coalescent mutation rate theta/2
        self.splitT = splitT
        self.mu = mu
        
        #PSMC parameters
        self.lh = lambdas#pairs of PSMC lambda_0 and lambda_1
        self.times = times
        self.numT = len(self.lh)#number of time intervals
        if len(self.times) != self.numT - 1:
            print("Unexpected number of time intervals")
            sys.exit(0)
        
        #Data parameters
        #Joint allele frequency spectrum: 0100,1100,0001,0101,1101,0011,0111
        self.dataJAFS = dataJAFS
        
        #Class variables
        self.lc = [[0,0] for i in range(self.numT)]#Corrected lambdas
        self.M = None#Differential equation matrix
        self.integralP = None#Integral of dif eq solution
        self.P0 = None#Initial condition for dif eq
        self.P1 = None#Values of solution at the end of the interval
        self.JAFS = [0 for i in range(7)]#Joint allele frequency spectrum: 0100,1100,0001,0101,1101,0011,0111
        
        #Class for EP size correction
        self.cl = CorrectLambda()
        print("MigrationInference class initialized.")
      
    def PrintError(self, func, text):
        func = func + "():"
        print("MigrationInference class error in function", func, text)
        sys.exit(0)
      
    def CorrectLambdas(self):
        p0 = [[1,0,0],[0,1,0]]
        if not self.correct:
            for t in range(self.numT):
                 self.lc[t][0],self.lc[t][1] = self.lh[t][0],self.lh[t][1]
            return True
        for t in range(self.splitT):
#            print(self.lh[t])
#            print(self.times[t])
#            print(p0)
            self.cl.SetInterval(self.lh[t], self.times[t], p0)
            try:
                sol = self.cl.SolveLambdaSystem()
            except optimize.nonlin.NoConvergence:
                return False
#            print("interval solution\t",sol)
            self.lc[t][0],self.lc[t][1] = sol[0][0],sol[0][1]
            p0 = sol[1]
        for t in range(self.splitT,self.numT):
            self.lc[t][0],self.lc[t][1] = (self.lh[t][0]+self.lh[t][1])/2,(self.lh[t][0]+self.lh[t][1])/2
        return True
            
    def JAFSpectrum(self):
        model = TwoPopulations(self.lc[0][0], self.lc[0][1], 1.0, 1.0)
        self.P0 = [0.0 for i in range( model.MSize() )]
        self.P0[2] = 1.0
        for interval in range(self.numT):
            print("Interval", interval)
            if interval < self.splitT:
                if interval == self.numT - 1 and self.mu[0] + self.mu[1] == 0:
                    self.PrintError("JAFSpectrum", "Infinite coalescent time. No migration.")
                model = TwoPopulations(self.lc[interval][0], self.lc[interval][1], self.mu[0], self.mu[1])
            else:
                model = OnePopulation(self.lc[interval][0])
            if interval == self.splitT:
                self.CollapsePops()
            self.M = model.SetMatrix()
            self.P0 = model.SetInitialConditions(self.P0)
            self.SolveDifEq(interval)
            self.P0 = model.UpdateInitialConditions(self.P1)
            if interval < self.numT - 1:
                self.integralP = model.UpdateIntegral(self.integralP, self.times[interval])
            for i in range( model.StateNum() ):
                jaf = model.StateToJAF(i)
                self.JAFS = [x + y*self.integralP[i] for x,y in zip(self.JAFS, jaf)]
    
    def PrintMatrix(self):
        matSize = self.M.shape[0]
        for i in range(matSize):
            for j in range(matSize):
                el = format(self.M.item(i,j), '.10g')
             #   if int(el) == 0:
             #       el = '.'
                print( el, end = "\t" )
            print("")
    
    def CollapsePops(self):
        Pc = [0 for i in range(8)]
        Pc[0] = sum(self.P0[0:9])
        Pc[1] = sum(self.P0[9:15])
        Pc[2] = sum(self.P0[15:23])
        Pc[3] = sum(self.P0[23:29])
        Pc[4] = sum(self.P0[29:33])
        Pc[5] = sum(self.P0[33:37])
        Pc[6] = sum(self.P0[37:41])
        Pc[7] = sum(self.P0[41:44])
        self.P0 = Pc
    
    def SolveDifEq(self, interval):
        sizeM = self.M.shape[0]
        if interval < self.numT - 1:
            T = self.times[interval]
            MET = linalg.expm( dot(self.M,T) )
            self.P1 = dot(MET,self.P0)
        else:
            self.P1 = [0 for i in range( len(self.P0) )]
        MI = linalg.inv(self.M)
        self.integralP = [x - y for x, y in zip(self.P1, self.P0)]
        self.integralP = dot(MI,self.integralP)
    
    def ObjectiveFunction(self, mu):
        self.cl.SetMu(mu[0], mu[1])
        res = self.CorrectLambdas()
        if not res:
            return -10**(100)
        if 0:
            print("ObjectiveFunction(): corrected values of lambdas are", self.lc)
        self.JAFSpectrum()
        norm = sum(self.JAFS)
        if 1:
            print("----------",self.JAFS[0]/norm,self.JAFS[1]/norm,sep="\t\t")
            print(self.JAFS[2]/norm,self.JAFS[3]/norm,self.JAFS[4]/norm,sep="\t\t")
            print(self.JAFS[5]/norm,self.JAFS[6]/norm,"----------",sep="\t\t")
            n = 1+1/2+1/3
            print("singletons", (self.JAFS[0]+self.JAFS[2])/norm, 1/n)
            print("doubletons", (self.JAFS[1]+self.JAFS[3]+self.JAFS[5])/norm, 1/(2*n))
            print("tripletons", (self.JAFS[4]+self.JAFS[6])/norm, 1/(3*n))
#        return 0
#        return self.Likelihood()
        llh = 0
        for i in range(7):
            llh += self.dataJAFS[i]*log(self.JAFS[i])
        print(llh)
        return -llh
    
    def Solve(self):
        print("Start solving the problem.")
        self.cl = CorrectLambda()
#        res = optimize.minimize(self.ObjectiveFunction, [0,0], method=’L-BFGS-B’, bounds = ((0, None), (0, None)))
#        print(res)
        #print(self.ObjectiveFunction([self.mu[0], self.mu[1]]))
#        optimize ObjectiveFunction(mu0, mu1)
    
    def Test(self):
        self.ObjectiveFunction([self.mu[0], self.mu[1]])
        

lambdas, times, mu, splitT, theta = [], [], [], None, None

l = [1,1/0.4]
#l = [1]
lambdas = []
for la in l:
    lambdas.append([la,la])
#print(lambdas)
#times = [2*0.125]
#times = []
#lambdas = [[1,3],[0.1,0.5],[2,4],[1,0.5],[10,10]]
#times = [2,3,1,5]
#mu = [0.0, 0.0]#[0.00001,0.00001]
#mu = [0.00001,0.00001]
#splitT = 3

#times = [2*0.125]
lambdas = [[1.0, 2.5], [1.4285714285714286, 1.4285714285714286]]
mu = [100.0, 300.0]
mu = [0.0, 0.0]
times = [2.5]
splitT = 1

lambdas = [[1.0, 2.5], [1.4285714285714286, 1.6666666666666667], [2.0, 2.0]]
mu = [100.0, 30.0]
times = [0.25, 0.35]
splitT = 2
dataJAFS = [0 for i in range(7)]

inputData = [times, lambdas]

if len(sys.argv) < 4:
    print("./migration <PSMC input file 1> <PSMC input file 2> <JAF spectrum file>")
    sys.exit(0)
fpsmc1 = sys.argv[1]
fpsmc2 = sys.argv[2]
fjafs  = sys.argv[3]
inputData = migrationIO.ReadPSMC(fpsmc1, fpsmc2)
dataJAFS = migrationIO.ReadJAFS(fjafs)

'''theta = 0.000001
lambdas = [[1.0, 1.0], [2.5, 2.5]]
mu = [0.01, 0.003]
times = [0.25]
splitT = 0'''
print("Input parameters:")
print("\tlambda = ", lambdas, sep = "")
print("\ttimes  = ", times, sep = "")
print("\tmigrat = ", mu, sep = "")
print("\tsplitT = ", splitT, sep = "")
print("\ttheta  = ", theta, sep = "")
print("End of input.")
for i in range( len(inputData[0]) - 1 ):
    splitT = i
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, False)
    Migration.ObjectiveFunction([mu[0], mu[1]])

#./scrm 4 1 -t 100000.0 -r 100000.0 250000000 -l 100000 -eN 0.0125 0.4 > sim.txt
#./scrm 4 1 -t 100000.0 -r 100000.0 250000000 -l 100000 -eN 0.0125 0.4 > sim.txt