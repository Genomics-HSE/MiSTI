#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
from math import log
from CorrectLambda import CorrectLambda
from MigrationInference import MigrationInference

def ReadPSMCFile(fn, RD = -1):
    maxRD = -1
    Tk = []
    Lk = []
    th = 0
    
    with open(fn) as f:
        for line in f:
            line = line.split()
            if line[0] == "RD":
                maxRD = int( line[1] )
        if maxRD == -1:
            print("Corrupted of empty input file")
            sys.exit(0)
        if RD == -1 or RD > maxRD:
            RD = maxRD

    with open(fn) as f:  
        for line in f:
            line = line.split()
            if line[0] != "RD" or int( line[1] ) != RD:
                continue
            while line[0] != "RS":
                if line[0] == "TR":
                    th = float(line[1])
                line = next(f)
                line = line.split()
            while line[0] != "PA":
                if line[0] != "RS":
                    print("Unexpected line.")
                    sys.exit(0)
                Tk.append( float(line[2]) )
                Lk.append( float(line[3]) )
                line = next(f)
                line = line.split()
            break
    data = [Tk, Lk, RD, th]
    return( data )

def ReadPSMC(fn1, fn2, RD = -1, doPlot = False):
    d1 = ReadPSMCFile(fn1, RD)
    d2 = ReadPSMCFile(fn2, RD)
    if d1[2] != d2[2]:
        print("Different RDs for input files 1 and 2.")
        sys.exit(0)
    d2[0] = [v*d2[3]/d1[3] for v in d2[0]]#rescale   time       by th1/th2
    d2[1] = [v*d2[3]/d1[3] for v in d2[1]]#rescale   epsize     by th2/th1 (compare with previous line!)
    Tk = []
    Lk1 = []
    Lk2 = []
    Tk = sorted( d1[0] + d2[0][1:] )
    
    j = 0
    for i in range( len(d1[0]) - 1 ):
        while Tk[j] < d1[0][i + 1]:
            Lk1.append( 1.0/d1[1][i] )
            j += 1
    while len(Lk1) < len(Tk):
        Lk1.append(1.0/d1[1][-1])
        
    j = 0
    for i in range( len(d2[0]) - 1 ):
        while Tk[j] < d2[0][i + 1]:
            Lk2.append( 1.0/d2[1][i] )
            j += 1
    while len(Lk2) < len(Tk):
        Lk2.append(1.0/d2[1][-1])
    scale = 1
    scale1 = 1
    if doPlot:
#    print("Here ready to plot")
        mu = 6.83e-8
        binsize = 100
        scale = d1[3]/(2.0*binsize*mu)
        scale1 = scale/2.0/1.0e4
        x = [v*scale for v in Tk]
        y1 = [scale1/v for v in Lk1]
        y2 = [scale1/v for v in Lk2]
#        plt.semilogx(x, y1, x, y2)
        AddToPlot(x, y1)
        AddToPlot(x, y2)
    L1tmp = [Lk1[i]]
    L2tmp = [Lk1[i]]
    Ttmp = [Tk[0]]
    '''    if True:
        for i in range( 1, len(Lk1) ):
            if Lk1[i] != Lk1[i-1] or Lk2[i] != Lk2[i-1]:
                L1tmp.append(Lk1[i])
                L2tmp.append(Lk1[i])
                Ttmp.append()'''
                
#        plt.savefig("temp2.png")
#    print("Here we are")
#    sys.exit(0)
    Lk = [[u, v] for u, v in zip(Lk1, Lk2)]
    Tk = [ u - v for u, v in zip(Tk[1:], Tk[:-1])]
    return( [Tk, Lk, scale, scale1] )
    
#    print(len(Tk))
#    print(len(Lk1))
#    for i in range(len(Tk)):
#        print(1/Lk1[i], "\t", Tk[i])

def OutputMigration(fout, mu, Migration):
    Migration.JAFSLikelyhood( mu )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times)+1)]
    
    outData = "#Migration ver 0.1\n"
    outData += "ST\t" + str(Migration.splitT) + "\n"#split times
    outData += "MU\t" + str(mu[0]) + "\t" + str(mu[1]) + "\n"#migration
    for i in range( len(times) ):
        outData += str(times[i]) + "\t" + str(Migration.lc[i][0]) + "\t" + str(Migration.lc[i][1]) + "\n"
    
    if fout == "":
        print(outData)
    else:
        fw = open(fout, 'w')
        fw.write(outData)
        fw.close()
    
def ReadMigration(fmigr, doPlot=False, scaleTime = 1, scaleEPS = 1):
    times = []
    lc1 = []
    lc2 = []
    splitT = None
    mu = None
    with open(fmigr) as f:
        line = next(f).rstrip()
        line = line.split(" ")
        print("Format version: ", line[2])
        line = next(f).rstrip()
        line = line.split("\t")
        splitT = int(line[1])
        line = next(f).rstrip()
        line = line.split("\t")
        mu = [float(line[1]), float(line[2])]
        for line in f:
            line = line.split("\t")
            times.append( float(line[0]) )
            lc1.append( float(line[1]) )
            lc2.append( float(line[2]) )
    if doPlot:
        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc1])
        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc2])
        splT=times[splitT]#sum(inputData[0][0:splitT])
        plt.axvline(splT*scaleTime, color='r')

def ReadJAFS(fn):
    jafs = []
    with open(fn) as f:
        for line in f:
            jafs.append( float(line) )
    if len(jafs) != 7:
        print("Unexpected number of lines in the JAFS file.")
        sys.exit(0)
    return(jafs)
    
def PlotInit(id=1):
    plt.figure(id)
    plt.semilogx()
    
def AddToPlot(times, lambdas, id=1):
    plt.figure(id)
    plt.step(times, lambdas)
    
def SavePlot(fout, id=1):
    plt.figure(id)
    plt.savefig(fout)