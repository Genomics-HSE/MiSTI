#!/usr/bin/env python3

#    Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)
#
#    This file is part of MiSTI.
#
#    MiSTI is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MiSTI is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MiSTI.  If not, see <https://www.gnu.org/licenses/>.


from scipy import (optimize,linalg,integrate)
from numpy import (dot,identity)
import numpy
import sys
from math import (exp,log,sqrt)



class CorrectLambda:
    def __init__(self):
        self.mixtureTH = 0.02
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
    
    def SetMixtureTH(self, th):
        self.mixtureTH = th
    
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
    
    def ExpectedCoalTimeOnePop(self, lam):
        return 1.0/lam-self.T/(exp(lam*self.T)-1)
        
    def ExpectedCoalTimeOnePopTmp(self, lam):
        pnc = exp(-lam*self.T)
        Tc = 1.0/lam-self.T/(1.0/pnc-1.0)
        return [Tc, pnc]
    
    def ExpectedCoalTimeOnePopNonConditional(self, lam):
        return (1-exp(-lam*self.T)*(1+lam*self.T))/lam
    
    def EPSFromExpectedCoalTime(self, Te, x0 = 1, prec = 1e-10):
        upperLimit = numpy.inf#10*self.lh[0]
        lowerLimit = 0.01*min(self.lh[0], self.lh[1])#0
        x1 = optimize.least_squares(lambda lam: self.ExpectedCoalTimeOnePop(lam) - Te, x0, bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
        return x1.x
        
    def FitSinglePop(self):
        pnc = [sum(self.P0[0]),sum(self.P0[1])]
        pnc = [pnc[0]/sum(pnc), pnc[1]/sum(pnc)]
        Te = pnc[0]*self.ExpectedCoalTimeOnePop(self.lh[0]) + pnc[1]*self.ExpectedCoalTimeOnePop(self.lh[1])
        return self.EPSFromExpectedCoalTime(Te, pnc[0]*self.lh[0] + pnc[1]*self.lh[1])
    
    def ExpectedCoalTimeTwoPop(self):
        Pnormed = [None, None]
        coalT = [[None, None], [None, None]]
        for i in [0,1]:
            Pnormed[i] = [v/sum(self.P0[i]) for v in self.P0[i]]
        Minv = linalg.inv(self.M)
        for i in [0,1]:
            mat = self.MET-identity(3)
            vec1 = dot(mat, Pnormed[i])
            vec1 = dot(Minv, dot(Minv, vec1))
            vec2 = dot(self.MET, Pnormed[i])
            pnc = sum(vec2)
            vec2 = dot(self.T, dot(Minv, vec2))
            vec = vec2 - vec1
            coalT[i][0] = (self.l[0]*vec[0] + self.l[1]*vec[1])/(1-pnc)
            coalT[i][1] = pnc
        return coalT
    
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
    
    def VariantRatio(self, Tc, pnc):
        singletons = pnc*2*self.T+(1.0-pnc)*2*Tc
        doubletons = (1.0-pnc)*(self.T-Tc)
        return singletons/doubletons
    
    def LambdaSystem(self,l):
        self.l = [l[0],l[1]]
        self.SetMatrix()
        self.MatrixExponent()
        coalT = self.ExpectedCoalTimeTwoPop()
        coalT1 = [self.ExpectedCoalTimeOnePopTmp(self.lh[0]), self.ExpectedCoalTimeOnePopTmp(self.lh[1])]
        return( coalT[0][0]-coalT1[0][0], coalT[1][0]-coalT1[1][0])
        
    def LambdaSystemRatio(self,l):#fit singletons to doubletons ratio
        self.l = [l[0],l[1]]
        self.SetMatrix()
        self.MatrixExponent()
        coalT = self.ExpectedCoalTimeTwoPop()
        coalT1 = [self.ExpectedCoalTimeOnePopTmp(self.lh[0]), self.ExpectedCoalTimeOnePopTmp(self.lh[1])]
        varRatio1pop = [self.VariantRatio(coalT1[0][0], coalT1[0][1]), self.VariantRatio(coalT1[1][0], coalT1[1][1])]
        varRatio2pop = [self.VariantRatio(coalT[0][0], coalT[0][1]), self.VariantRatio(coalT[1][0], coalT[1][1])]
        return( varRatio1pop[0]-varRatio2pop[0], varRatio1pop[1]-varRatio2pop[1])
 
    def LambdaSystem1(self,l):
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

    def SolveNoMigration1(self):
        lc = [None, None]
        A1 = self.P0[0][0]/sum(self.P0[0])
        A2 = self.P0[0][1]/sum(self.P0[0])
        A3 = self.P0[1][0]/sum(self.P0[1])
        A4 = self.P0[1][1]/sum(self.P0[1])
        C1 = self.P0[0][2]/sum(self.P0[0])
        C2 = self.P0[1][2]/sum(self.P0[1])
        D  =  A1*A4-A2*A3
        B1 =  A4/D
        B2 = -A2/D
        B3 = -A3/D
        B4 =  A1/D
        X1 = exp(-self.lh[0]*self.T)-C1
        X2 = exp(-self.lh[1]*self.T)-C2
        if B1*X1+B2*X2 > 0 and B3*X1+B4*X2 > 0:
            lc[0] = -log(B1*X1+B2*X2)/self.T
            lc[1] = -log(B3*X1+B4*X2)/self.T
        else:
            lc = [-1, -1]
        p0 = [self.P0[0][0]*exp(-lc[0]*self.T), self.P0[0][1]*exp(-lc[1]*self.T), self.P0[0][2]]
        p1 = [self.P0[1][0]*exp(-lc[0]*self.T), self.P0[1][1]*exp(-lc[1]*self.T), self.P0[1][2]]
        return [lc,[p0,p1]]
    
    def LambdaSystemNoMigration(self,l):
        coalT = [None, None]
        for i in [0,1]:
            pnc = self.pr0[i][0]*exp(-l[0]*self.T)+self.pr0[i][1]*exp(-l[1]*self.T)+self.pr0[i][2]
            #            print("pnc[",i,"]=",pnc)
            coalT[i] = (self.pr0[i][0]*self.ExpectedCoalTimeOnePopNonConditional(l[0])+self.pr0[i][1]*self.ExpectedCoalTimeOnePopNonConditional(l[1]))/(1-pnc)
            '''        print("l=",l)
        print("lh=",self.lh)
        print("time int=",self.T)
        print("coalT=",coalT)
        print("exptT=",[self.ExpectedCoalTimeOnePop(self.lh[0]),self.ExpectedCoalTimeOnePop(self.lh[1])])
        print(self.pr0[0])
        print(self.pr0[1])
        sys.exit(0)'''
        return( coalT[0]- self.ExpectedCoalTimeOnePop(self.lh[0]), coalT[1]- self.ExpectedCoalTimeOnePop(self.lh[1]) )
    
    def SolveNoMigration(self):
        self.pr0 = [None, None]
        for i in [0,1]:
            self.pr0[i] = [v/sum(self.P0[i]) for v in self.P0[i]]
        upperLimit = numpy.inf#10*self.lh[0]
        lowerLimit = 0.01*min(self.lh[0], self.lh[1])#0
        prec = 1e-10
        x1 = optimize.least_squares(self.LambdaSystemNoMigration, self.lh, bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
        lc = x1.x
        p0 = [self.P0[0][0]*exp(-lc[0]*self.T), self.P0[0][1]*exp(-lc[1]*self.T), self.P0[0][2]]
        p1 = [self.P0[1][0]*exp(-lc[0]*self.T), self.P0[1][1]*exp(-lc[1]*self.T), self.P0[1][2]]
        return [lc,[p0,p1]]

    def SolveLambdaSystem(self, cpfit = True, prec = 1e-10, normEps=0.02):
        mixture = 0
        for i in range(3):
            mixture += (self.P0[0][i]/sum(self.P0[0])-self.P0[1][i]/sum(self.P0[1]))**2
        mixture = sqrt(mixture)
        if mixture < self.mixtureTH:
            return [[-1, -1], self.P0]
        if self.mu[0] + self.mu[1] < prec:
            if not cpfit:
                return self.SolveNoMigration()
            else:
                return self.SolveNoMigration1()
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
        stretch = True
        if stretch:
            timeTmp = self.T
            self.T = self.T/timeTmp
            self.mu = [self.mu[0]*timeTmp, self.mu[1]*timeTmp]
            self.lh = [self.lh[0]*timeTmp, self.lh[1]*timeTmp]
        upperLimit = numpy.inf#10*self.lh[0]
        lowerLimit = 0.001*min(self.lh[0], self.lh[1])#0
        if not cpfit:
            x1 = optimize.least_squares(self.LambdaSystem, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
        else:
            x1 = optimize.least_squares(self.LambdaSystem1, [self.lh[0],self.lh[1]], bounds = (lowerLimit, upperLimit), gtol = prec, xtol = prec)
        x = x1.x
        self.l = x
        if stretch:
            self.T = timeTmp
            self.mu = [self.mu[0]/timeTmp, self.mu[1]/timeTmp]
            self.lh = [self.lh[0]/timeTmp, self.lh[1]/timeTmp]
            self.l =  [self.l[0]/timeTmp,  self.l[1]/timeTmp]
        self.SetMatrix()
        self.MatrixExponent()
        p0 = dot(self.MET,self.P0[0])
        p1 = dot(self.MET,self.P0[1])
        return [self.l,[p0,p1]]
        
#if __name__ == "__main__":
    
    
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
