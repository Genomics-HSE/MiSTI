#!/usr/bin/env python3

#Copyright (c) 2018 Vladimir Shchur (vlshchur@gmail.com)

import sys
import collections
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)

from CorrectLambda import CorrectLambda

#lineage = collections.namedtuple('lineage', ['d0', 'd1', 'pop' ])#d0: number of descendents in population 0; d1: number of descendents in population 1; pop: current population of lineage

class lineage:
    def __init__(self, d0, d1):
        self.d0 = d0
        self.d1 = d1

class OnePopulation:
    def __init__(self, l1):
        self.la = l1
        self.Msize = 8
    
    def MSize(self):
        return( self.Msize )
    
    def StateNum(self):
        return( self.Msize )
    
    def PrintError(self, func, text):
        func = func + "():"
        print("TwoPopulations class error in function", func, text)
        sys.exit(0)
        
    def CheckState(self, state):
        tn0, tn1 = 0, 0
        state[:] = sorted(state, key=lambda line: line.d0, reverse = True)
        state[:] = sorted(state, key=lambda line: line.d0+line.d1, reverse=True)
        for val in state:
            tn0 += val.d0
            tn1 += val.d1
        if tn0 != 2 or tn1 != 2:
            print("CheckState() not passed: expected number of lineages is 2 and 2, instead there are", tn0, "and", tn1, "lineages")
            sys.exit(0)
    
    def MapStateToInd(self, state):
        self.CheckState(state)
        if len(state) == 4:
            return 0
        if len(state) == 3:
            if state[0].d0 == 2:
                return 1
            elif state[0].d1 == 1:
                return 2
            elif state[0].d0 == 0:
                return 3
        if len(state) == 2:
            if state[0].d0 == 2 and state[0].d1 == 1:
                return 4
            elif state[0].d0 == 1 and state[0].d1 == 2:
                return 5
            elif state[0].d0 == 2 and state[0].d1 == 0:
                return 6
            elif state[0].d0 == 1 and state[0].d1 == 1:
                return 7
        return self.Msize
    
    def MapIndToState(self, ind):
        if not isinstance(ind, int) or ind < 0 or ind >= self.Msize:
            print("Unexpected index value", ind, ", index should be an integer between 0 and", self.Msize - 1, ".")
            sys.exit(0)
        if ind == 0:
            state = [lineage(1, 0), lineage(1, 0), lineage(0, 1), lineage(0, 1)]
        elif ind == 1:
            state = [lineage(2, 0), lineage(0, 1), lineage(0, 1)]
        elif ind == 2:
            state = [lineage(1, 1), lineage(1, 0), lineage(0, 1)]
        elif ind == 3:
            state = [lineage(0, 2), lineage(1, 0), lineage(1, 0)]
        elif ind == 4:
            state = [lineage(2, 1), lineage(0, 1)]
        elif ind == 5:
            state = [lineage(1, 2), lineage(1, 0)]
        elif ind == 6:
            state = [lineage(2, 0), lineage(0, 2)]
        else:
            state = [lineage(1, 1), lineage(1, 1)]
        self.CheckState(state)
        return state
    
    def StateToJAF(self, sti):
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
    
    def PrintState(self, state):
        st = ""
        for val in state:
            st = st + "(" + str(val.d0) + "," + str(val.d1) + ") "
        return st
    
    def SetInitialConditions( self , P0 ):
        return( P0 )
    
    def UpdateInitialConditions( self , P0):
        return( P0 )
        
    def UpdateIntegral(self, integralP, T):
        return( integralP )
    
    def SetMatrix(self):
        MM = [[0 for col in range(self.Msize)] for row in range(self.Msize)]#transition matrix
        for i in range(self.Msize):
            self.UpdateMatrixCol(i, MM)
        MM = mat(MM)
        return(MM)
    
    def UpdateMatrixCol(self, ind, MM):
        state = self.MapIndToState(ind)
        total = 0.0
#        self.la = 1.0
        for i in range(len(state)):
            #Coalescence
            for j in range(i+1, len(state)):
                if j == i:
                    continue
                state1 = list(state)
                for k in sorted([i,j], reverse=True):
                    del state1[k]
                line = lineage(state[i].d0+state[j].d0, state[i].d1+state[j].d1)
                state1.append(line)
                ind2 = self.MapStateToInd(state1)
                if ind2 != self.Msize:
                    MM[ind2][ind] += self.la
                total += self.la
        MM[ind][ind] -= total
        
    def Test(self):
        map_test_passed = True
        for i in range(self.Msize):
            st = self.MapIndToState(i)
            if 0:
                print(i, self.PrintState(st), self.MapStateToInd(st), sep="\t")
            if self.MapStateToInd(st) != i:
                map_test_passed = False
        if map_test_passed:
            print("Map test passed succesfully, map is identity.")
        else:
            self.PrintError("Test", "map test failed.")
            
        for i in range(self.Msize):
            st = self.MapIndToState(i)