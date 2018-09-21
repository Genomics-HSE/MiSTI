#!/Users/shchur/anaconda3/bin/python3

import sys
import os
import collections
from scipy import (linalg,optimize)
from numpy import (dot,identity,mat)
from math import exp
import multiprocessing
import matplotlib.pyplot as plt
import random
import time

import argparse

def func(var):
    with open("african/Din_San.jaf") as f:
        test = "text"
    print(test)

class TestClass:
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.l = [1, 2, 3]

M = 3e9
m1 = 342071
m2 = 374671
m3 = 107940
m4 = 98035
m5 = 97138
m6 = 85478
m7 = 72114
T1 = (m1/2+m3-(2*m6+m5)*(6*m7+m5)/8/m5)/M
T2 = (m2/2+m4-(2*m7+m5)*(6*m6+m5)/8/m5)/M

mu = 1.25e-8
Ygen = 30
print(T1/mu*Ygen)
print(T2/mu*Ygen)

if 0:
    timesMS = [0, 0.0075, 0.0225, 0.125, 0.5, 1.875]
    epsMS = [[0.2, 0.2], [0.2, 0.45], [0.5, 0.5], [0.3, 0.3], [0.5, 0.5], [1.3, 1.3]]
    timesTmp = [0]
    epsTmp = []
    for i in range(len(timesMS)-1):
        maxTime = timesTmp[-1]
        deltaT = timesMS[i+1] - timesMS[i]
        timesTmp += [maxTime + (j+1)*deltaT/10.0 for j in range(10)]
        epsTmp += [epsMS[i] for j in range(10)]
    epsTmp.append(epsMS[-1])
    timesMS = timesTmp
    epsMS = epsTmp
    #for u,v in zip(timesMS, epsMS):
    #    print(u, "\t", v)
    times = [2*(u-v) for u, v in zip(timesMS[1:], timesMS[:-1])]
    eps = [[1.0/u[0], 1.0/u[1]] for u in epsMS]



sys.exit(0)
