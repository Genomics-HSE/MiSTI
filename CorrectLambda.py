#!/Users/shchur/anaconda3/bin/python3
from scipy import (optimize,linalg)
from numpy import (dot,identity)
from math import (exp,log)

class CorrectLambda:
#    def __init__(self):
#        self.lh = [lh0, lh1]#PSMC inferred lambda
#        self.mu = [mu0, mu1]
#        self.P0 = [[1,0,0],[0,1,0]]
#        self.T = T
    
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
    
    def LambdaEqaution(self, npop):
        assert npop == 0 or npop == 1, "Population number should be 0 or 1."
        if self.T == -1:
            nch = sum(self.P0[npop])
        else:
            nch = (1-exp(-self.lh[npop]*self.T))*sum(self.P0[npop])
        self.ComputeExpectation(npop)
        nc = self.l[0]*self.Pexp[0]+self.l[1]*self.Pexp[1]
        return nc-nch
    
    def LambdaSystem(self,l):
        self.l = [l[0],l[1]]
        self.SetMatrix()
        self.MatrixExponent()
        return([self.LambdaEqaution(0), self.LambdaEqaution(1)])
#        print(self.LambdaEqaution(0))
#        print(self.LambdaEqaution(1))
    
    def SolveLambdaSystem(self, prec = 1e-14):
#        x = optimize.broyden1(self.LambdaSystem, [0.00001,0.00001], f_tol=prec)
        x = optimize.broyden1(self.LambdaSystem, [self.lh[0],self.lh[1]], f_tol=prec)
        self.l = x
        self.SetMatrix()
        self.MatrixExponent()
        p0 = dot(self.MET,self.P0[0])
        p1 = dot(self.MET,self.P0[1])
        return [x,[p0,p1]]