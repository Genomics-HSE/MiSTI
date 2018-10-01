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

class lineage:#d0: number of descendents in population 0; d1: number of descendents in population 1; pop: current population of lineage
    def __init__(self, d0, d1, pop):
        self.d0 = d0
        self.d1 = d1
        self.pop = pop

class TwoPopulations:
    def __init__(self, l1, l2, m1, m2):
        if m1 < 0 or m2 < 0:
            self.PrintError("TwoPopulations", "migration rates cannot be negative")
        if l1 < 0 or l2 < 0:
            self.PrintError("TwoPopulations", "coalescent rates cannot be negative")
        self.mu = [m1, m2]
        self.la = [l1, l2]
        self.P0 = []
        self.stationary = []
        self.Msize = 44
        if m1 + m2 == 0:
            for i in range(self.Msize):
                st = self.MapIndToState(i)
                if len(st) == 2 and st[0].pop != st[1].pop:
                    self.stationary.append(i)
            if False:
                print("Stationary states:", self.stationary)
        self.Test()
    
    def MSize(self):
        return self.Msize - len(self.stationary)
    
    def StateNum(self):
        return(self.Msize)
    
    def PrintError(self, func, text):
        func = func + "():"
        print("TwoPopulations class error in function", func, text)
        sys.exit(0)
        
    def CheckState(self, state):
        tn0, tn1 = 0, 0
        state[:] = sorted(state, key=lambda line: line.pop)
        state[:] = sorted(state, key=lambda line: line.d0, reverse = True)
        state[:] = sorted(state, key=lambda line: line.d0+line.d1, reverse=True)
        for val in state:
            tn0 += val.d0
            tn1 += val.d1
        if tn0 != 2 or tn1 != 2:
            text = "CheckState() not passed: expected number of lineages is 2 and 2, instead there are " + str(tn0) + " and " + str(tn1) + " lineages"
            self.PrintError("CheckState", text)
    
    def MapStateToInd(self, state):
        self.CheckState(state)
        if len(state) == 4:
            i, j = 0, 0
            for val in state:
                if val.d0 == 0:
                    i += val.pop
                else:
                    j += val.pop
            return i+3*j
        if len(state) == 3:
            if state[0].d0 + state[0].d1 != 2:
                text = "Unexpected number of lineages, expected 2, recieved " + str(state[0].d0 + state[0].d1)
                self.PrintError("MapStateToInd", text)
            if state[0].d0 == 2:
                return 9 + 3*state[0].pop + state[1].pop + state[2].pop
            elif state[0].d0 == 1:
                return 9 + 6 + 4*state[0].pop + 2*state[1].pop + state[2].pop
            elif state[0].d0 == 0:
                return 9 + 6 + 8 + 3*state[0].pop + state[1].pop + state[2].pop
        if len(state) == 2:
            if state[0].d0 == 2 and state[0].d1 == 1:
                return 29 + 2*state[0].pop + state[1].pop
            elif state[0].d0 == 1 and state[0].d1 == 2:
                return 29 + 4 + 2*state[0].pop + state[1].pop
            elif state[0].d0 == 2 and state[0].d1 == 0:
                return 29 + 4 + 4 + 2*state[0].pop + state[1].pop
            elif state[0].d0 == 1 and state[0].d1 == 1:
                return 29 + 4 + 4 + 4 + state[0].pop + state[1].pop
        return self.Msize
    
    def MapIndToState(self, ind):
        if not isinstance(ind, int) or ind < 0 or ind >= self.Msize:
            text = "Unexpected index value " + str(ind) + ", index should be an integer between 0 and " + str(self.Msize - 1) + "."
            self.PrintError("MapIndToState", text)
        if ind < 9:
            state = [lineage(1, 0, 0), lineage(1, 0, 0), lineage(0, 1, 0), lineage(0, 1, 0)]
            i = ind//3
            j = ind%3
            while i > 0:
                self.UpdateLineagePop(state, i-1, 1)
                i -= 1
            while j > 0:
                self.UpdateLineagePop(state, 2+j-1, 1)
                j -= 1
        elif ind < 29:
            if ind < 15:
                state = [lineage(2, 0, 0), lineage(0, 1, 0), lineage(0, 1, 0)]
                self.UpdateLineagePop(state, 0, (ind - 9)//3)
                i = (ind - 9)%3
                while i > 0:
                    self.UpdateLineagePop(state, i, 1)
                    i -= 1
            elif ind < 23:
                state = [lineage(1, 1, 0), lineage(1, 0, 0), lineage(0, 1, 0)]
                self.UpdateLineagePop(state, 0, (ind - 15)//4)
                self.UpdateLineagePop(state, 1, ((ind - 15)%4)//2 )
                self.UpdateLineagePop(state, 2, ((ind - 15)%4)%2 )
            else:
                state = [lineage(0, 2, 0), lineage(1, 0, 0), lineage(1, 0, 0)]
                self.UpdateLineagePop(state, 0, (ind - 23)//3)
                i = (ind - 23)%3
                while i > 0:
                    self.UpdateLineagePop(state, i, 1)
                    i -= 1
        else:
            if ind < 33:
                state = [lineage(2, 1, 0), lineage(0, 1, 0)]
                popTmp = ind - 29
                self.UpdateLineagePop(state, 0, popTmp//2)
                self.UpdateLineagePop(state, 1, popTmp%2)
            elif ind < 37:
                state = [lineage(1, 2, 0), lineage(1, 0, 0)]
                popTmp = ind - 33
                self.UpdateLineagePop(state, 0, popTmp//2)
                self.UpdateLineagePop(state, 1, popTmp%2)
            elif ind < 41:
                state = [lineage(2, 0, 0), lineage(0, 2, 0)]
                self.UpdateLineagePop(state, 0, (ind - 37)//2)
                self.UpdateLineagePop(state, 1, (ind - 37)%2)
            else:
                state = [lineage(1, 1, 0), lineage(1, 1, 0)]
                if ind == 43:
                    self.UpdateLineagePop(state, 0, 1)
                if ind == 42 or ind == 43:
                    self.UpdateLineagePop(state, 1, 1)
        self.CheckState(state)
        return state
    
    def StateToJAF(self, sti):
        #sti += sum(i <= sti for i in self.stationary)
#        for el in self.stationary:
#            if sti >= el:
#                sti += 1
        state = self.MapIndToState(sti)
        jaf = [0 for i in range(7)]
        for lineage in state:
            if lineage.d0 == 0:
                if lineage.d1 == 1:
                    jaf[2] += 1
                elif lineage.d1 == 2:
                    jaf[5] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
            if lineage.d0 == 1:
                if lineage.d1 == 0:
                    jaf[0] += 1
                elif lineage.d1 == 1:
                    jaf[3] += 1
                elif lineage.d1 == 2:
                    jaf[6] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
            if lineage.d0 == 2:
                if lineage.d1 == 0:
                    jaf[1] += 1
                elif lineage.d1 == 1:
                    jaf[4] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
        return(jaf)
    
    def UpdateLineagePop(self, state, ind, pop):
        l = lineage(state[ind].d0, state[ind].d1, pop)
        state[ind] = l
    
    def PrintState(self, state):
        st = ""
        for val in state:
            st = st + "(" + str(val.d0) + "," + str(val.d1) + "," + str(val.pop) + ") "
        return st
        
    def SetMatrix(self):
        MM = [[0 for col in range(self.Msize)] for row in range(self.Msize)]#transition matrix
        for i in range(self.Msize):
            self.UpdateMatrixCol(i, MM)
        MM = numpy.delete(MM, self.stationary, 0)
        MM = numpy.delete(MM, self.stationary, 1)
        MM = mat(MM)
        return(MM)
        
    def SetInitialConditions( self , P0 ):
        self.P0 = P0
        if self.mu[0] + self.mu[1] == 0:
            P0 = numpy.delete(P0, self.stationary)
        return( P0 )
    
    def AncientSampleP0(self, P0):
        newP0 = [0.0 for v in P0]
        for i in range(self.Msize):
            st = self.MapIndToState(i)
            check = 0
            for state in st:
                if state.d0 == 1 and state.d1 == 0 and state.pop == 0:
                    check += 1
            if check == 2:
                newP0[2] += P0[i]
            check = 0
            for state in st:
                if state.d0 == 2 and state.d1 == 0 and state.pop == 0:
                    check += 1
            if check == 1:
                newP0[11] += P0[i]
        return(newP0)
    
    def UpdateInitialConditions( self , P0):
        if len(P0) == self.Msize:
            return( P0 )
        if len(P0) != self.Msize - len(self.stationary):
            text = "unexpected length of initial conditions vector " + len(P0)
            self.PrintError("UpdateInitialConditions", text)
        for ind in self.stationary:
            P0 = numpy.insert( P0, ind, 0 )
        nP0 = [v for v in P0]
        for ind in self.stationary:
            st = self.MapIndToState(ind)
            c = [st[0].d0*st[0].pop + st[1].d0*st[1].pop, st[0].d1*st[0].pop + st[1].d1*st[1].pop]
            for i in range(self.Msize):
                st1 = self.MapIndToState(i)
                c1 = [0,0]
                for lineage in st1:
                    c1[0] += lineage.d0*lineage.pop
                    c1[1] += lineage.d1*lineage.pop
                if c1[0] == c[0] and c1[1] == c[1]:
                    nP0[ind] += self.P0[i] - P0[i]
        for i in range(self.Msize):
            if nP0[i] > 0:
                st = self.MapIndToState(i)
        return( nP0 )
    
    def UpdateIntegral(self, integralP, T):
        if len(integralP) == self.Msize:
            return( integralP )
        if len(integralP) != self.Msize - len(self.stationary):
            text = "unexpected length of initial conditions vector " + len(P0)
            self.PrintError("UpdateIntegral", text)
        for ind in self.stationary:
            integralP = numpy.insert( integralP, ind, 0 )
        niP = [v for v in integralP]
        for ind in self.stationary:
            st = self.MapIndToState(ind)
            c = [st[0].d0*st[0].pop + st[1].d0*st[1].pop, st[0].d1*st[0].pop + st[1].d1*st[1].pop]
            for i in range(self.Msize):
                st1 = self.MapIndToState(i)
                c1 = [0,0]
                for lineage in st1:
                    c1[0] += lineage.d0*lineage.pop
                    c1[1] += lineage.d1*lineage.pop
                if c1[0] == c[0] and c1[1] == c[1]:
                    niP[ind] += T*self.P0[i] - integralP[i]
        return( niP )
    
    def PrintMatrix(self):
        matSize = self.Msize - len(self.stationary)
        for i in range(matSize):
            for j in range(matSize):
                el = format(self.MM.item(i,j), '.10g')
                if float(el) == 0:
                    el = '.'
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
                pr = self.PrintState(st)
                if int(pr) == 0:
                    pr = '.'
                print ( pr )
    
    def UpdateMatrixCol(self, ind, MM):
        state = self.MapIndToState(ind)
        total = 0.0
        for i in range(len(state)):
            #Migration
            state1 = list(state)
            self.UpdateLineagePop(state1, i, 1-state[i].pop )
            ind2 = self.MapStateToInd(state1)
            MM[ind2][ind] += self.mu[ state[i].pop ]
            total += self.mu[ state[i].pop ]
            #Coalescence
            for j in range(i+1, len(state)):
                if state[j].pop != state[i].pop:
                    continue
                state1 = list(state)
                for k in sorted([i,j], reverse=True):
                    del state1[k]
                line = lineage(state[i].d0+state[j].d0, state[i].d1+state[j].d1, state[i].pop)
                state1.append(line)
                ind2 = self.MapStateToInd(state1)
                if ind2 != self.Msize:
                    MM[ind2][ind] += self.la[ state[i].pop ]
                total += self.la[ state[i].pop ]
        MM[ind][ind] -= total
    
    def Test(self):
        return
        map_test_passed = True
        if 0:
            for i in range(self.Msize):
                st = self.MapIndToState(i)
                print(self.PrintState(st))
                if self.MapStateToInd(st) != i:
                    map_test_passed = False
            return
        if map_test_passed:
            print("Map test passed succesfully, map is identity.")
        else:
            print("Map test failed.")
        if 0:
           jaf = self.StateToJAF(i)
        #               j = '' + `int(jaf[0])` + `int(jaf[1])` + `int(jaf[2])` + `int(jaf[3])`
           print(self.PrintState(st), jaf)
        if 1:
           self.mu = [1, 10]
           self.la = [0, 0]
           self.MM = self.SetMatrix()
           self.PrintMatrix()
           sys.exit(0)
           
#tp = TwoPopulations(0, 0, 1, 100)