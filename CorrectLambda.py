#!/usr/bin/env python3

#Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)

from scipy import (optimize,linalg,integrate)
from numpy import (dot,identity)
import numpy
import sys
from math import (exp,log,sqrt)



class CorrectLambda:
#    def __init__(self):
#        self.doBroyden1 = True
#        self.lh = [lh0, lh1]#PSMC inferred lambda
#        self.mu = [mu0, mu1]
#        self.P0 = [[1,0,0],[0,1,0]]
#        self.T = T
    COUNT = 0
    
    def PrintError(self, func, text):
        func = func + "():"
        print("MigrationInference class error in function", func, text)
        sys.exit(0)
    
    def SetMu(self, mu0, mu1):
        self.mu = [mu0, mu1]
        
    def SetInterval(self, lh, T, P0):#T=-1 means T=infinity
        self.lh = [lh[0], lh[1]]
        self.T = T
        self.P0 = P0
    
    def SetMatrix(self):
        self.M = [[-2*self.mu[0]-self.l[0], 0, self.mu[1]],[0, -2*self.mu[1]-self.l[1], self.mu[0]],[2*self.mu[0], 2*self.mu[1], -self.mu[0] - self.mu[1]]]
        
    def MatrixExponent(self):
        if self.T == -1:
            self.MET =[[0,0,0],[0,0,0],[0,0,0]]
        else:
            self.MET = linalg.expm( dot(self.M,self.T) )#exponent of M times T
        
    def ComputeExpectation(self, npop):
        self.Pexp = dot(linalg.inv(self.M), dot(self.MET-identity(3),self.P0[npop]))
    
    def LambdaEquation1(self, npop):
        assert npop == 0 or npop == 1, "Population number should be 0 or 1."
        if self.T == -1:
            self.PrintError("LambdaEquation", "lambda correction for the last interval is not implemented")
            nch = sum(self.P0[npop])
        else:
            nch = (1-exp(-self.lh[npop]*self.T))*sum(self.P0[npop])
        self.ComputeExpectation(npop)
        nc = self.l[0]*self.Pexp[0]+self.l[1]*self.Pexp[1]
        return nc-nch
        
    def LambdaEquation(self, npop):
        assert npop == 0 or npop == 1, "Population number should be 0 or 1."
        if self.T == -1:
            self.PrintError("LambdaEquation", "lambda correction for the last interval is not implemented")
            nch = sum(self.P0[npop])
        else:
            nch = exp(-self.lh[npop]*self.T)*sum(self.P0[npop])
        p0 = dot(self.MET,self.P0[npop])
        nc = sum(p0)
        return nc-nch
    
    def LambdaSystem(self,l):
        self.l = [l[0],l[1]]
        self.SetMatrix()
        self.MatrixExponent()
        return( numpy.array([self.LambdaEquation(0), self.LambdaEquation(1)]) )
#        print(self.LambdaEqaution(0))
#        print(self.LambdaEqaution(1))

    def CoalRateInterval(self,l):
        self.l = [l[0],l[1]]
        self.SetMatrix()
        print(self.M)
        self.ODE([1, 0, 0],0)
    
    def ODE(self, y, t):
        return( dot(self.M, y) )
    
    def CoalRates(self, intervals, splitT, discr = 100):#interval = [time, lambda1, lambda2, mu1, mu2], discr = number of intervals in the discretization
        p0 = [1.0, 0.0, 0.0]
        solution = [[],[]]
        for i in range( len(intervals) ):
            inter = intervals[i]
            t1 = inter[0]
            if inter[0] < splitT:
                self.SetMu(inter[3], inter[4])
                self.l = [1.0/inter[1],1.0/inter[2]]
                self.SetMatrix()
                t2 = intervals[i+1][0]
                times = [t1+(t2-t1)/discr*p for p in range(discr + 1)]
                sol = integrate.odeint(self.ODE, p0, times)
                p0 = sol[-1]
                for i in range( len(sol) ):
                    prob = list(sol[i])
                    rate = (self.l[0]*prob[0] + self.l[1]*prob[1])/sum(prob)
                    solution[0].append(rate)
                    solution[1].append(times[i])
        return(solution)

    def StationarySystem(self, l):
        eq0 = -(2*self.mu[0]+l[0]) + self.mu[1]/self.a[0]
        eq1 = -(2*self.mu[1]+l[1]) + self.mu[0]/self.a[1]
        eq2 = 2*self.mu[0]*self.a[0] + 2*self.mu[1]*self.a[1] - (self.mu[0]+self.mu[1])
        return([eq0-eq2, eq1-eq2])

    def SolveLambdaSystemExperimental(self, prec = 1e-10, normEps=0.02):#computes lambdas from stationary distribution when |self.P0[0]-self.P0[1]| < eps. The issue: how to capture changes in ef. pop. sizes?
        normV0 = 0
        normV1 = 0
        normD = 0
        for i in range(3):
            normV0 += ( self.P0[0][i] )**2
            normV1 += ( self.P0[1][i] )**2
            normD += (self.P0[0][i] - self.P0[1][i])**2
        normV0 = sqrt(normV0)
        normV1 = sqrt(normV1)
        normD = sqrt(normD)
        x = None
        if normD < normEps*min(normV0, normV1):
            statV = [(el1+el2)/2.0 for el1, el2 in zip(self.P0[0], self.P0[1])]
            upperLimit = numpy.inf
            lowerLimit = 0
            self.a = [statV[0]/statV[2], statV[1]/statV[2]]
            x1 = optimize.least_squares(self.StationarySystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
            x = x1.x
        else:
            upperLimit = numpy.inf#10*self.lh[0]
            lowerLimit = 0#0.1*self.lh[0]
            x1 = optimize.least_squares(self.LambdaSystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
            x = x1.x
        self.l = x
        self.SetMatrix()
        self.MatrixExponent()
        p0 = dot(self.MET,self.P0[0])
        p1 = dot(self.MET,self.P0[1])
        return [x,[p0,p1]]

    def SolveLambdaSystem(self, prec = 1e-10, normEps=0.02):
        normV0 = 0
        normV1 = 0
        normD = 0
        for i in range(3):
            normV0 += ( self.P0[0][i] )**2
            normV1 += ( self.P0[1][i] )**2
            normD += (self.P0[0][i] - self.P0[1][i])**2
        normV0 = sqrt(normV0)
        normV1 = sqrt(normV1)
        normD = sqrt(normD)
        if normD < normEps*min(normV0, normV1):
            nlh = (self.lh[0]+self.lh[1])/2.0
            self.lh[0],self.lh[1] = nlh,nlh
        x = None
#        x = optimize.broyden1(self.LambdaSystem, [self.lh[0],self.lh[1]], f_tol=prec)
        upperLimit = numpy.inf#10*self.lh[0]
        lowerLimit = 0.1*min(self.lh[0], self.lh[1])#0
        x1 = optimize.least_squares(self.LambdaSystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
        x = x1.x
        '''        if False:
            try:
                x = optimize.broyden1(self.LambdaSystem, [self.lh[0],self.lh[1]], f_tol=prec)
            except:
                self.doBroyden1 = False
                upperLimit = numpy.inf#10*self.lh[0]
                lowerLimit = 0.1*self.lh[0]#0
                x1 = optimize.least_squares(self.LambdaSystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
                x = x1.x
        if False:
            upperLimit = numpy.inf#10*self.lh[0]
            lowerLimit = 0.1*self.lh[0]#0
            x1 = optimize.least_squares(self.LambdaSystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
            x = x1.x'''
        self.l = x
        self.SetMatrix()
        self.MatrixExponent()
        p0 = dot(self.MET,self.P0[0])
        p1 = dot(self.MET,self.P0[1])
        '''        print("\tmatrix=\n",self.MET)
        print(self.P0[0])
        print(self.P0[1])
        print("\tlc=", x)
        print("\tlh=", self.lh)
        print("\tmu=", self.mu)
        print("\tt=", self.T)
        print("\t\tLambdaSystemSolution=", self.LambdaSystem(x))'''
        return [x,[p0,p1]]
        
        
'''
cl = CorrectLambda()
inter1 = [0, 1, 1, 1, 1]
inter2 = [0.02, 0.05, 0.05, 1, 1]
inter3 = [0.075, 0.5, 0.5, 1, 1]
inter4 = [2.5, 1, 1, 1, 1]
splitT = 2.5
intervals = [inter1, inter2, inter3, inter4]
coalRates = cl.CoalRates(intervals, splitT)
print(coalRates)'''