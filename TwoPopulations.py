#!/Users/shchur/anaconda3/bin/python3

import sys
import collections
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
import math
from math import (exp,log)

#lineage = collections.namedtuple('lineage', ['d0', 'd1', 'pop' ])#d0: number of descendents in population 0; d1: number of descendents in population 1; pop: current population of lineage

class lineage:
    def __init__(self, d0, d1, pop):
        self.d0 = d0
        self.d1 = d1
        self.pop = pop

class TwoPopulations:
    def __init__(self, l1, l2, m1, m2):
        self.mu = [m1, m2]
        self.la = [l1, l2]
        self.Msize = 44
        
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
            print("CheckState() not passed: expected number of lineages is 2 and 2, instead there are", tn0, "and", tn1, "lineages")
            sys.exit(0)
    
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
                print("Unexpected number of lineages, expected 2, recieved", state[0].d0 + state[0].d1)
                sys.exit(0)
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
            print("Unexpected index value", ind, ", index should be an integer between 0 and", self.Msize - 1, ".")
            sys.exit(0)
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
        MM = mat(MM)
        return(MM)
    
    def UpdateMatrixCol(self, ind, MM):
        state = self.MapIndToState(ind)
        total = 0.0
        for i in range(len(state)):
            #Migration
            state1 = list(state)
            self.UpdateLineagePop(state1, i, (state[i].pop+1)%2)
            ind2 = self.MapStateToInd(state1)
            MM[ind][ind2] += self.mu[ state[i].pop ]
            total += self.mu[ state[i].pop ]
            #Coalescence
            for j in range(i+1, len(state)):
                if j == i or state[j].pop != state[i].pop:
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
        map_test_passed = True
        for i in range(self.Msize):
            st = self.MapIndToState(i)
            if self.MapStateToInd(st) != i:
                map_test_passed = False
        if map_test_passed:
            print("Map test passed succesfully, map is identity.")
        else:
            print("Map test failed.")