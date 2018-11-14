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


import sys
import collections
import numpy
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)
import time
import multiprocessing

from CorrectLambda import CorrectLambda
from TwoPopulations import TwoPopulations
from OnePopulation import OnePopulation

class MigrationInference:
    COUNT_LLH = 0
    CORRECTION_CALLED = 0
    CORRECTION_FAILED = 0
    LLH_CONST = 0
    
    
    def __init__(self, times, lambdas, dataJAFS, mu, splitT, **kwargs):#thrh for theta and rho
        self.debug = False
        self.enableOutput = False
        if "debug" in kwargs:
            if kwargs["debug"]:
                self.debug = True
                self.enableOutput = True
        if "enableOutput" in kwargs:
            if kwargs["enableOutput"]:
                self.enableOutput = True
        if self.enableOutput:
            print("MigrationInference: output enabled.")
        self.cpfit = False
        if "cpfit" in kwargs:
            if kwargs["cpfit"]:
                self.cpfit = True
        
        self.correct = True
        if "trueEPS" in kwargs:
            if kwargs["trueEPS"]:
                self.correct = False
        self.smooth = False
        if "smooth" in kwargs:
            if kwargs["smooth"]:
                self.smooth = True

        self.unfolded = False
        if "unfolded" in kwargs:
            if kwargs["unfolded"]:
                self.unfolded = True
        
        self.thrh = [1.0, 1.0]#theta and rho from PSMC1
        if "thrh" in kwargs:
            if len(kwargs["thrh"]) == 2:
                self.thrh = kwargs["thrh"]
        
        self.sampleDate = 0#Dating of the second sample, by default it is 0.0 - present time
        if "sampleDate" in kwargs:
            self.sampleDate = kwargs["sampleDate"]#time in generation units
            if self.sampleDate > 0:
                print("Sample date needs to be fixed in this version.")
                sys.exit(0)
        
        if splitT < self.sampleDate:
            self.PrintError("__init__", "cannot initialise class with split time being more recent than sample date.")
        
        #Model parameters
        splitFraction = splitT%1
        splitT = int(splitT)
        if splitT - 1 > len(times):
            self.PrintError("__init__", "Invalid value for split time, cannot create Migration class instance.")
        if splitFraction != 0.0:
            t1 = splitFraction*times[splitT]
            t2 = times[splitT] - t1
            times[splitT] = t1
            times.insert(splitT + 1, t2)
            lambdas.insert(splitT + 1, lambdas[splitT])
            splitT += 1
        
        #PSMC parameters
        self.lh = list(lambdas)#pairs of PSMC lambda_0 and lambda_1
        self.times = times
        self.numT = len(self.lh)#number of time intervals
        if len(self.times) != self.numT - 1:
            print("Unexpected number of time intervals")
            sys.exit(0)
        '''discr_step = 0.0001 
        self.discr = 200
        timeMap = [0]
        if discr_step > 0:
            timesTmp = []
            lhTmp = []
            discr = self.discr
            newST = self.discr*splitT
            splitT = 120
            for i in range( splitT ):
                discr = math.ceil(self.times[i]/discr_step)
                discr = self.discr
                timeMap.append(discr+timeMap[-1])
                timesTmp += [self.times[i]/discr for j in range(discr)]
                lhTmp += [self.lh[i] for j in range(discr)]
            timeMap += list(range(len(timesTmp),len(timesTmp)+len(self.times[splitT:])))
            #newST = len(timesTmp)
            lhTmp += list(self.lh[splitT:])
            timesTmp += list(self.times[splitT:])
            #lhTmp.append(self.lh[-1])
            self.times = timesTmp
            self.lh = lhTmp
            self.numT = len(self.lh)
            splitT = newST
            for i in range(len(mu)):
                mu[i][1] = timeMap[int(mu[i][1])]
                mu[i][2] = timeMap[int(mu[i][2])]'''

        self.discr = 100
        if self.discr > 1:
            timesTmp = []
            lhTmp = []
            discr = self.discr
            splitT = splitT*discr
            if False:
                for i in range(len(self.times)):
                    self.times[i] = self.times[i]*discr
                for i in range(len(self.lh)):
                    self.lh[i] = [self.lh[i][0]/discr,self.lh[i][1]/discr]
            for i in range( self.numT - 1 ):
                timesTmp += [self.times[i]/discr for j in range(discr)]
                lhTmp += [self.lh[i] for j in range(discr)]
            lhTmp.append(self.lh[-1])
            self.times = timesTmp
            self.lh = lhTmp
            self.numT = len(self.lh)

            for i in range(len(mu)):
                mu[i][1] = int(mu[i][1])*discr
                mu[i][2] = int(mu[i][2])*discr
                if False:
                    mu[i][3] = float(mu[i][3])/discr

        self.splitT = splitT
        self.mu = list(self.lh)#initialize migration with the same size as lambdas
        self.SetModel(mu)
        
        #Data parameters
        #Joint allele frequency spectrum: 0100,1100,0001,0101,1101,0011,0111
        self.snps = dataJAFS.jafs[0]
        self.dataJAFS = [el for el in dataJAFS.jafs[1:]]
        
        if MigrationInference.LLH_CONST == 0:
            if self.unfolded:
                if dataJAFS.ufLLHConst == None:
                    for j in range(1, self.snps+1):
                        MigrationInference.LLH_CONST += log(j)
                    for i in range(7):
                        for j in range(1, self.dataJAFS[i]+1):
                            MigrationInference.LLH_CONST -= log(j)
                else:
                    MigrationInference.LLH_CONST = dataJAFS.ufLLHConst
            else:
                if dataJAFS.fLLHConst == None:
                    for j in range(1, self.snps+1):
                        MigrationInference.LLH_CONST += log(j)
                    for j in range(1, self.dataJAFS[0]+self.dataJAFS[6]+1):
                        MigrationInference.LLH_CONST -= log(j)
                    for j in range(1, self.dataJAFS[1]+self.dataJAFS[5]+1):
                        MigrationInference.LLH_CONST -= log(j)
                    for j in range(1, self.dataJAFS[2]+self.dataJAFS[4]+1):
                        MigrationInference.LLH_CONST -= log(j)
                    for j in range(1, self.dataJAFS[3]+1):
                        MigrationInference.LLH_CONST -= log(j)
                else:
                    MigrationInference.LLH_CONST = dataJAFS.fLLHConst
        
        #Class variables
        self.lc = [[1,1] for i in range(self.numT)]#Corrected lambdas
        self.M = None#Differential equation matrix
        self.integralP = None#Integral of dif eq solution
        self.P0 = None#Initial condition for dif eq
        self.P1 = None#Values of solution at the end of the interval
        self.JAFS = [0 for i in range(7)]#Joint allele frequency spectrum: 0100,1100,0001,0101,1101,0011,0111
        
        #Class for EP size correction
        self.cl = CorrectLambda()
        if "mixtureTH" in kwargs:
            self.cl.SetMixtureTH(kwargs["mixtureTH"])
#        self.cl.SetMu(mu[0], mu[1])
        
        #Plotting options - TODO
        self.doPlot = False
        if "doPlot" in kwargs:
            if kwargs["doPlot"]:
                self.doPlot = True
                self.scaleEPS = 1
                self.scaleT = 1
                if "scaleEPS" in kwargs:
                    self.scaleEPS = float(kwargs[scaleEPS])
                if "scaleT" in kwargs:
                    self.scaleT = float(kwargs[scaleT])
        if self.debug:
            print("MigrationInference class initialized. Class size", self.numT)
    
    def SetModel(self, params):
        self.optPars = []
        for i in range(len(self.mu)):
            self.mu[i] = [0.0,0.0]
        for el in params:
            popInd = int(el[0]) - 1
            if popInd != 0 and popInd != 1:
                self.PrintError("SetModel", "Population index should be 1 or 2.")
            migStart = int(el[1])
            migEnd = int(el[2])
            migVal = float(el[3])
            migVar = int(el[4])
            for i in range(migStart, migEnd):
                if self.mu[i][popInd] != 0:
                    self.PrintError("SetModel", "Migration rate intervals should not overlap.")
                self.mu[i][popInd] = migVal
            if migVar == 1:
                self.optPars.append([popInd, migStart, migEnd, migVal])
        self.optParsSize = len(self.optPars)
    
    def MapParameters(self, params):
        if len(params) != self.optParsSize:
            self.PrintError("MapParameters", "Incorrect number of parameters.")
        for i in range(self.optParsSize):
            for j in range(self.optPars[i][1], self.optPars[i][2]):
                self.mu[j][self.optPars[i][0]] = params[i]
    
    def PrintError(self, func, text):
        func = func + "():"
        print("MigrationInference class error in function", func, text)
        sys.exit(0)
      
    def CorrectLambdas(self):
        MigrationInference.CORRECTION_CALLED += 1
        p0 = [[1,0,0],[0,1,0]]
        self.Pr = [[[1.0,0.0],[0.0,1.0],[0.0,0.0]]]
        nc = [0, 0]#Probability for not coalescing
        for t in range(self.splitT):
#            print(self.lh[t])
#            print(self.times[t])
#            print(p0)
            self.cl.SetMu(self.mu[t][0], self.mu[t][1])
            if not self.correct:# or self.mu[t][0] + self.mu[t][1] == 0:
                self.lc[t][0],self.lc[t][1] = self.lh[t][0],self.lh[t][1]
            else:
                self.cl.SetInterval(self.lh[t], self.times[t], p0)
                try:
                    sol = self.cl.SolveLambdaSystem(self.cpfit)
                except optimize.nonlin.NoConvergence:
                    print("lh=", self.lh[t])
                    print("t=", self.times[t])
                    print("mu=", self.mu[t])
                    print("p0=", p0)
                    MigrationInference.CORRECTION_FAILED += 1
                    sys.exit(0)
                    return False
                if self.enableOutput:
                    print("lh=", self.lh[t])
                    print("lc=", sol[0])
                    print("t=", self.times[t])
                    print("initial conditions", p0[0], "\t", p0[1])
#                print("time = ", t, "\tinterval solution\t",sol)
                self.lc[t][0],self.lc[t][1] = sol[0][0],sol[0][1]
                if sol[0][0] < 0 or sol[0][1] < 0:
                    MigrationInference.CORRECTION_FAILED += 1
                    return False
                p0 = sol[1]
            self.Pr.append([[p0[0][0],p0[1][0]],[p0[0][1],p0[1][1]],[p0[0][2],p0[1][2]]])
#            nc[0] += -self.times[t]*self.lh[t][0]
#            nc[1] += -self.times[t]*self.lh[t][1]
            nc[0] = sum(p0[0])
            nc[1] = sum(p0[1])
        for t in range(self.splitT,self.numT - 1):
#            self.lc[t][0],self.lc[t][1] = (self.lh[t][0]+self.lh[t][1])/2,(self.lh[t][0]+self.lh[t][1])/2
            if self.times[t] == 0:
                self.lc[t][0], self.lc[t][1] = 1, 1
                continue
            if not self.cpfit:
                self.cl.SetInterval(self.lh[t], self.times[t], [[exp(nc[0]), 0, 0], [exp(nc[1]), 0, 0]])
                lam = self.cl.FitSinglePop()[0]
                self.lc[t][0] = lam
                self.lc[t][1] = lam
            else:
                pnc = ( exp(-self.times[t]*self.lh[t][0]) + exp(nc[1] - nc[0] - self.times[t]*self.lh[t][1]) )/( 1 + exp( nc[1] - nc[0] ) )
                self.lc[t][0] = -log(pnc)/self.times[t]
                self.lc[t][1] = -log(pnc)/self.times[t]
            nc[0] += -self.times[t]*self.lc[t][0]
            nc[1] += -self.times[t]*self.lc[t][1]
        t = self.numT - 1
        pr0 = exp(nc[0])
        pr1 = exp(nc[1])
        lam = (pr0+pr1)/(pr0/self.lh[t][0]+pr1/self.lh[t][1])
        self.lc[t][0] = lam
        self.lc[t][1] = lam
        self.Smooth()
        return True
    
    def Smooth(self):
        if not self.smooth:
            return
#        self.SmoothSplitTime(2, 2)
        self.SmoothConst(0)
        self.SmoothConst(1)
    
    def SmoothConst(self, indiv):
        k = 0
        lam = self.lh[0][indiv]
        time = 0.0
        nc = 0.0
        while k < self.splitT:
            j = k
            while abs(self.lh[j][indiv]-lam) < 1e-10 and j < self.numT-1:
                nc += self.lc[j][indiv]*self.times[j]
                time += self.times[j]
                j += 1
                if j == self.splitT:
                    break
            for i in range(k, j):
                self.lc[i][indiv] = nc/time
            lam = self.lh[j][indiv]
            nc = 0.0
            time = 0.0
            k = j

    def SmoothInterval(self, interval, left, right = None):
        if not self.smooth:
            return
        if right == None:
            right = left
        if interval - left < 0:
            text = "Unexpected value for left smoothing: " + str(left) + ", value cannot be larger than " + str(interval)
            self.PrintError("SmoothInterval", text)
        if interval + right - 1 >= self.numT:
            text = "Unexpected value for right smoothing: " + str(right) + ", value cannot be larger than " + str(self.numT - interval + 1)
            self.PrintError("SmoothInterval", text)
        ncSmInt = [0, 0]
        sumT = 0
        for k in [0, 1]:
            for i in range(interval - left, interval + right):
                ncSmInt[k] += -self.lc[i][k]*self.times[i]
        for i in range(interval - left, interval + right):
            sumT += self.times[i]
        nc = [0, 0]
        for k in [0, 1]:
            for i in range(self.splitT - left):
                nc[k] += -self.lc[i][k]*self.times[i]
        pnc = (1 + exp(nc[1] - nc[0] + ncSmInt[1] - ncSmInt[0])) / (1 + exp(nc[1] - nc[0]))
        pnc = ncSmInt[0] + log(pnc)
        lsmoothed = - pnc / sum(times[self.splitT - left:self.splitT + right])
        lsmoothed = - ncSmInt[0]/sumT
        for k in [0, 1]:
            for i in range(self.splitT - left, self.splitT + right):
                self.lc[i][k] = lsmoothed
    
    def SmoothSplitTime(self, left, right = None):
        if not self.smooth:
            return
        if right == None:
            right = left
        if self.splitT - left < 0:
            text = "Unexpected value for left smoothing: " + str(left) + ", value cannot be larger than " + str(self.splitT)
            self.PrintError("SmoothSplitTime", text)
        if self.splitT + right - 1 >= self.numT:
            text = "Unexpected value for right smoothing: " + str(right) + ", value cannot be larger than " + str(self.numT - self.splitT + 1)
            self.PrintError("SmoothSplitTime", text)
        ncSmInt = [0, 0]
        sumT = 0
        for k in [0, 1]:
            for i in range(self.splitT - left, self.splitT + right):
                ncSmInt[k] += -self.lc[i][k]*self.times[i]
        for i in range(self.splitT - left, self.splitT + right):
            sumT += self.times[i]
        nc = [0, 0]
        for k in [0, 1]:
            for i in range(self.splitT - left):
                nc[k] += -self.lc[i][k]*self.times[i]
        pnc = (1 + exp(nc[1] - nc[0] + ncSmInt[1] - ncSmInt[0])) / (1 + exp(nc[1] - nc[0]))
        pnc = ncSmInt[0] + log(pnc)
        lsmoothed = - pnc / sum(times[self.splitT - left:self.splitT + right])
        lsmoothed = - ncSmInt[0]/sumT
        for k in [0, 1]:
            for i in range(self.splitT - left, self.splitT + right):
                self.lc[i][k] = lsmoothed
    
    def JAFSpectrum(self):
        model = TwoPopulations(self.lc[0][0], self.lc[0][1], 1.0, 1.0)
        self.P0 = [0.0 for i in range( model.MSize() )]
        self.P0[2] = 1.0
        pnc = 1#used from present time to ancient genome time
        for interval in range(self.numT):
            if interval < self.splitT:
                if interval == self.numT - 1 and self.mu[interval][0] + self.mu[interval][1] == 0.0:
                    self.PrintError("JAFSpectrum", "Infinite coalescent time. No migration.")
                model = TwoPopulations(self.lc[interval][0], self.lc[interval][1], self.mu[interval][0], self.mu[interval][1])
            else:
                model = OnePopulation(self.lc[interval][0])
            if interval == self.sampleDate:
                model_tmp = TwoPopulations(1, 1, 0, 0)
                self.P0 = model_tmp.AncientSampleP0(self.P0)
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
                if interval < self.sampleDate:
                    for i in range(2, len(jaf)):
                        jaf[i] = 0
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
        
    def JAFSLikelyhood(self, mu):
        MigrationInference.COUNT_LLH += 1
        self.llh = -10**9
        for v in mu:
            if v < 0:
                return self.llh#float('-inf')
#        self.mu[0],self.mu[1]=mu[0],mu[1]
        self.MapParameters(mu)
        res = self.CorrectLambdas()
        if not res:
            return self.llh#float('-inf') # -10**(10)
        if self.enableOutput:
            print("JAFSLikelyhood():   initial values of lambdas are ", self.lh)
            print("JAFSLikelyhood(): corrected values of lambdas are ", self.lc)
        self.JAFSpectrum()
        norm = sum(self.JAFS)
        self.JAFS = [v/norm for v in self.JAFS]
        if self.debug:
            print("----------",self.JAFS[0],self.JAFS[1],sep="\t\t")
            print(self.JAFS[2],self.JAFS[3],self.JAFS[4],sep="\t\t")
            print(self.JAFS[5],self.JAFS[6],"----------",sep="\t\t")
            
            print("----------",self.dataJAFS[0]/sum(self.dataJAFS),self.dataJAFS[1]/sum(self.dataJAFS),sep="\t\t")
            print(self.dataJAFS[2]/sum(self.dataJAFS),self.dataJAFS[3]/sum(self.dataJAFS),self.dataJAFS[4]/sum(self.dataJAFS),sep="\t\t")
            print(self.dataJAFS[5]/sum(self.dataJAFS),self.dataJAFS[6]/sum(self.dataJAFS),"----------",sep="\t\t")
            n = 1+1/2+1/3
            print("singletons", (self.JAFS[0]+self.JAFS[2]), 1/n)
            print("doubletons", (self.JAFS[1]+self.JAFS[3]+self.JAFS[5]), 1/(2*n))
            print("tripletons", (self.JAFS[4]+self.JAFS[6]), 1/(3*n))
            print("JAFS = ", self.JAFS)
#        return 0
#        return self.Likelihood()
        llh = MigrationInference.LLH_CONST
        if not self.unfolded:
            llh += (self.dataJAFS[0]+self.dataJAFS[6])*log(self.JAFS[0]+self.JAFS[6])
            llh += (self.dataJAFS[1]+self.dataJAFS[5])*log(self.JAFS[1]+self.JAFS[5])
            llh += (self.dataJAFS[2]+self.dataJAFS[4])*log(self.JAFS[2]+self.JAFS[4])
            llh += self.dataJAFS[3]*log(self.JAFS[3])
        else:
            for i in range(7):
#            print("self.dataJAFS[i]", self.dataJAFS[i], "\t\tlog(self.JAFS[i])", log(self.JAFS[i]), "\t\tself.JAFS[i]", self.JAFS[i])
                llh += self.dataJAFS[i]*log(self.JAFS[i])
#        print("full log llh=", llh)
#        for i in range(self.numT):
#            print(i, "\t", self.lc[i][0], "\t", self.lc[i][1])
        self.llh = llh
        return( llh )
    
    def ObjectiveFunction(self, mu):
        return( -self.JAFSLikelyhood( mu ) )
    
    def Solve(self, tol=1e-4, globalOpt = False):
        if self.optParsSize > 0:
            mu0 = [val[3] for val in self.optPars]
            if globalOpt:
                res = optimize.basinhopping(self.ObjectiveFunction, mu0, T=0.5, minimizer_kwargs=dict(method='Nelder-Mead'))
            else:
                res = optimize.minimize(self.ObjectiveFunction, mu0, method='Nelder-Mead', options={'xatol': tol, 'fatol': tol, 'maxiter': 100, 'disp': True })
        #res = optimize.minimize(self.ObjectiveFunction, mu0, method='BFGS', options={'gtol': tol })
            return([res.x, -res.fun])
        else:
            return([[], self.JAFSLikelyhood([])])
        
    def Report():
#        print("Split time ", self.splitT)
        print("Total number of likelihood function calls is", MigrationInference.COUNT_LLH)
        print("Lambda correction called", MigrationInference.CORRECTION_CALLED, "times.")
        print("Lambda correction failed", MigrationInference.CORRECTION_FAILED, "times.")
