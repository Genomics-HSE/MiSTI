#!/Users/shchur/anaconda3/bin/python3

import sys
import collections
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)

from CorrectLambda import CorrectLambda
from TwoPopulations import TwoPopulations
from OnePopulation import OnePopulation

class MigrationInference:
    def __init__(self, lambdas, times, mu, splitT, theta, correct = True):
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
        
        #Class variables
        self.lc = [[0,0] for i in range(self.numT)]#Corrected lambdas
        self.M = None#Differential equation matrix
        self.integralP = None#Integral of dif eq solution
        self.P0 = None#Initial condition for dif eq
        self.P1 = None#Values of solution at the end of the interval
        self.JAFS = [0 for i in range(7)]#Joint allele frequency spectrum: 0100,1100,0001,0101,1101,0011,0111
    
    def PrintMatrix(self):
        for i in range(self.Msize):
            for j in range(self.Msize):
                el = format(self.MM.item(i,j), '.10g')
                print( el, end = "\t" )
            print("")
            
    def PrintMatrixRow(self, rn):
        st = self.MapIndToState(rn)
        print (self.PrintState(st) )
        for i in range(self.Msize):
            el = format(self.MM.item(rn, i), '.10g')
            print( el, end = "\t" )
        print("")
        for i in range(self.Msize):
            if self.MM.item(rn, i) != 0:
                st = self.MapIndToState(i)
                print (self.PrintState(st) )
        
    def CorrectLambdas(self):
        p0 = [[1,0,0],[0,1,0]]
        if not self.correct:
            for t in range(self.numT):
                 self.lc[t][0],self.lc[t][1] = self.lh[t][0],self.lh[t][1]
            return
        for t in range(self.splitT):
#            print(self.lh[t])
#            print(self.times[t])
#            print(p0)
            self.cl.SetInterval(self.lh[t], self.times[t], p0)
            try:
                sol = self.cl.SolveLambdaSystem()
            except optimize.nonlin.NoConvergence:
                return False
            print("interval solution\t",sol)
            self.lc[t][0],self.lc[t][1] = sol[0][0],sol[0][1]
            p0 = sol[1]
        for t in range(self.splitT,self.numT):
            self.lc[t][0],self.lc[t][1] = (self.lh[t][0]+self.lh[t][1])/2,(self.lh[t][0]+self.lh[t][1])/2
        return True
            
    def JAFSpectrum(self):
        model = TwoPopulations(self.lc[0][0], self.lc[0][1], self.mu[0], self.mu[1])
        self.P0 = [0.0 for i in range(model.Msize)]
        self.P0[2] = 1.0
        for interval in range(self.numT):
            print("Interval", interval)
            if interval < self.splitT:
                model = TwoPopulations(self.lc[interval][0], self.lc[interval][1], self.mu[0], self.mu[1])
            else:
                model = OnePopulation(self.lc[interval][0])
            if interval == self.splitT:
                self.CollapsePops()
            self.M = model.SetMatrix()
#            print(self.M)
            self.SolveDifEq(interval)
            for i in range(model.Msize):
                jaf = model.StateToJAF(i)
                self.JAFS = [x + y*self.integralP[i] for x,y in zip(self.JAFS, jaf)]
            self.P0 = self.P1
    
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
        print(self.P0)
        if interval < self.numT - 1:
            T = self.times[interval]
            MET = linalg.expm( dot(self.M,T) )
            self.P1 = dot(MET,self.P0)
        else:
            self.P1 = [0 for i in range( len(self.P0) )]
        MI = linalg.inv(self.M)
        sizeM = self.M.shape[0]
        self.integralP = [x - y for x, y in zip(self.P1, self.P0)]
        self.integralP = dot(MI,self.integralP)
    
    def ObjectiveFunction(self, mu):
        self.cl.SetMu(mu[0], mu[1])
        res = self.CorrectLambdas()
        if not res:
            return -10**(100)
        print(self.lc)
        self.JAFSpectrum()
        norm = sum(self.JAFS)
        print("----------",self.JAFS[0]/norm,self.JAFS[1]/norm,sep="\t\t")
        print(self.JAFS[2]/norm,self.JAFS[3]/norm,self.JAFS[4]/norm,sep="\t\t")
        print(self.JAFS[5]/norm,self.JAFS[6]/norm,"----------",sep="\t\t")
        n = 1+1/2+1/3
        print("singletons", (self.JAFS[0]+self.JAFS[2])/norm, 1/n)
        print("doubletons", (self.JAFS[1]+self.JAFS[3]+self.JAFS[5])/norm, 1/(2*n))
        print("tripletons", (self.JAFS[4]+self.JAFS[6])/norm, 1/(3*n))
#        return self.Likelihood()
        return log(1)
    
    def Solve(self):
        self.cl = CorrectLambda()
        self.ObjectiveFunction([self.mu[0], self.mu[1]])
        #print(self.ObjectiveFunction([self.mu[0], self.mu[1]]))
#        optimize ObjectiveFunction(mu0, mu1)
    
    def Test(self):
        self.Solve()
        

lambdas, times, mu, splitT, theta = [], [], [], None, None

l = [1,1/0.4]
#l = [1]
lambdas = []
for la in l:
    lambdas.append([la,la])
print(lambdas)
#times = [2*0.125]
#times = []
#lambdas = [[1,3],[0.1,0.5],[2,4],[1,0.5],[10,10]]
#times = [2,3,1,5]
#mu = [0.0, 0.0]#[0.00001,0.00001]
#mu = [0.00001,0.00001]
#splitT = 3

times = [2*0.125]
splitT = 3

theta = 0.000001
print(lambdas)
print(times)
print(mu)
print(splitT)
print(theta)
Migration = MigrationInference(lambdas, times, mu, splitT, theta, False)
Migration.Test()

#./scrm 4 1 -t 100000.0 -r 100000.0 250000000 -l 100000 -eN 0.0125 0.4 > sim.txt
#./scrm 4 1 -t 100000.0 -r 100000.0 250000000 -l 100000 -eN 0.0125 0.4 > sim.txt