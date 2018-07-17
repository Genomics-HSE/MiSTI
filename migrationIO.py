#!/usr/bin/env python3
import sys
import matplotlib.pyplot as plt
from math import log
from CorrectLambda import CorrectLambda
from MigrationInference import MigrationInference
import argparse

class MigData:
    def __init__():
        self.splitT = None
        self.migStart = None
        self.migEnd  = None
        self.times = None
        self.lambda1 = None
        self.lambda2 = None
        self.thrh = None
        if "splitT" in kwargs:
            self.splitT = kwargs["splitT"]
        if "migStart" in kwargs:
            self.migStart = kwargs["migStart"]
        if "migEnd" in kwargs:
            self.migEnd = kwargs["migEnd"]
        if "times" in kwargs:
            self.times = kwargs["times"]
        if "lambda1" in kwargs:
            self.lambda1 = kwargs["lambda1"]
        if "lambda2" in kwargs:
            self.lambda2 = kwargs["lambda2"]
        if "thrh" in kwargs:
            self.thrh = kwargs["thrh"]

def SetScaling():
    mu = 1.1e-8#6.83e-8
    binsize = 100
    scaling = [mu, binsize]
    return(scaling)

def PrintErr(*args, sep="", endl="\n"):
    message = ""
    for word in args:
        message += str(word) + sep
    message += endl
    sys.stderr.write(message)

def ReadPSMCFile(fn, RD = -1):
    maxRD = -1
    Tk = []
    Lk = []
    th = 0
    rh = 0
    
    with open(fn) as f:
        for line in f:
            line = line.split()
            if line[0] == "RD":
                maxRD = int( line[1] )
        if maxRD == -1:
            print("Corrupted or empty input file")
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
                    rh = float(line[2])
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
    data = [Tk, Lk, RD, th, rh]
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
    scaling = SetScaling()
    mu = scaling[0]
    binsize = scaling[1]
    scale = d1[3]/(2.0*binsize*mu)
    scale1 = scale/2.0/1.0e4
    if doPlot:
#    print("Here ready to plot")
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
    return( [Tk, Lk, scale, scale1, d1[3], d1[4]] )#time, coalescent rates, N_0 (assuming default bin size = 100), effective population size/10000 rescale factor, theta and rho (from PSMC)
    
#    print(len(Tk))
#    print(len(Lk1))
#    for i in range(len(Tk)):
#        print(1/Lk1[i], "\t", Tk[i])

def OutputMigration(fout, mu, Migration):
    llh = Migration.JAFSLikelyhood( mu )
    print("llh = ", llh)
#    print( vars(Migration) )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times)+1)]   
    outData = "#Migration ver 0.3\n"
    outData += "ST\t" + str(Migration.splitT) + "\n"#split times
    outData += "MS\t" + str(Migration.migStart) + "\n"#migration start
    outData += "ME\t" + str(Migration.migEnd) + "\n"#migration end
    outData += "MU\t" + str(mu[0]) + "\t" + str(mu[1]) + "\n"#migration
    outData += "TR\t" + str(Migration.thrh[0]) + "\t" + str(Migration.thrh[1]) + "\n"#migration
    outData += "SFS\t" + str(0) + "\t" + "\t".join(map(str, Migration.JAFS)) + "\n"#expected SFS
    for i in range( len(times) ):
        outData += "RS\t" + str(times[i]) + "\t" + str(Migration.lc[i][0]) + "\t" + str(Migration.lc[i][1]) + "\n"
    
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
    thrh = [1.0, 1.0]
    splitT = None
    mu = None
    migStart = None
    migEnd = None
    with open(fmigr) as f:
        line = next(f).rstrip()
        line = line.split(" ")
        version = float(line[2])
        print("Format version: ", version)
        if version < 0.3:
            PrintErr("File version is not supported anymore.")
            sys.exit(0)
        
        line = next(f).rstrip()
        line = line.split("\t")
        splitT = int(line[1])
        
        line = next(f).rstrip()
        line = line.split("\t")
        migStart = int(line[1])
        
        line = next(f).rstrip()
        line = line.split("\t")
        migEnd = int(line[1])
        
        line = next(f).rstrip()
        line = line.split("\t")
        mu = [float(line[1]), float(line[2])]
        
        line = next(f).rstrip()
        line = line.split("\t")
        thrh = [float(line[1]), float(line[2])]
        
        line = next(f).rstrip()
        line = line.split("\t")
        jaf = [float(line[1]), float(line[2]), float(line[3]), float(line[4]), float(line[5]), float(line[6]), float(line[7])]
        
        for line in f:
            line = line.split("\t")
            times.append( float(line[0])*scaleTime )
            lc1.append( float(line[1])/scaleEPS )
            lc2.append( float(line[2])/scaleEPS )        
    if doPlot:
        lc1 = [1.0/max(v*scaleEPS,0.1)*scaleEPS for v in lc1]
        lc2 = [1.0/max(v*scaleEPS,0.1)*scaleEPS for v in lc2]
#        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc1])
#        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc2])
        AddToPlot(times, lc1)
        AddToPlot(times, lc2)
        splT=times[splitT]#sum(inputData[0][0:splitT])
        plt.axvline(splT, color='k', alpha=0.1)
    data = [splitT, migStart, migEnd, times, lc1, lc2, thrh]
    return(data)

def ReadJAFS(fn):
    jafs = []
    with open(fn) as f:
        line = next(f).rstrip()
        line = line.split(" ")
        print("Format version: ", line[2])
        line = next(f).rstrip()
        line = line.split("\t")
        jafs.append( int(line[1]) )
        for line in f:
            line = line.split("\t")
            jafs.append( int(line[1]) )
    if len(jafs) != 8:
        print("Unexpected number of lines in the JAFS file.")
        sys.exit(0)
    return(jafs)

def ReadMS(argument_string):
    PrintErr("WARNING: ReadMS() is not safe to use, the function has many assumptions on the ms command line")
    args = argument_string.split(" ")
    genMS = [[[0.0, 1.0]], [[0.0, 1.0]]]
    splitT = 0
    for i in range( len(args) ):
        if args[i] == "-n" and False:#FIXME
            pop = int(args[i+1]) - 1
            genMS[pop].append([0.0, float(args[i+2])])
        if args[i] == "-en":
            pop = int(args[i+2]) - 1
            genMS[pop].append([float(args[i+1]), float(args[i+3])])
        if args[i] == "-eN":
            genMS[0].append([float(args[i+1]), float(args[i+2])])
            genMS[1].append([float(args[i+1]), float(args[i+2])])
        if args[i] == "-ej":
            splitT = float( args[i+1] )
    splitT = splitT*2.0
    scale = 1
    scale1 = 1
    theta = 1#FIXME
    scaling = SetScaling()
    mu = scaling[0]
    binsize = scaling[1]
    scale = theta/(2.0*binsize*mu)
    scale1 = scale/2.0/1.0e4
    Tk, Lk1, Lk2 = [], [], []
    for el in sorted(genMS[0], key=lambda val: val[0]):
        Tk.append(el[0]*2.0)
        Lk1.append(1.0/el[1])
    for el in sorted(genMS[1], key=lambda val: val[0]):
        Lk2.append(1.0/el[1])
    splitT = [i for i,x in enumerate(Tk) if x == splitT]
    Lk = [[u, v] for u, v in zip(Lk1, Lk2)]
    Tk = [ u - v for u, v in zip(Tk[1:], Tk[:-1])]
    return( [Tk, Lk, scale, scale1, splitT[0]] )
#4 1000 -t 8196 -r 1355 3000000 -l -I 2 2 2 -n 2 1.0 -em 0.0 1 2 2.0 -em 0.0 2 1 2.0 -en 0.01 1 0.05 -en 0.01 2 0.05 -en 0.0375 1 0.5 -en 0.0375 2 0.5 -ej 1.25 2 1 -eM 1.25 0.0 -eN 1.25 1.0

def PlotInit(id=1):
    plt.figure(id)
    plt.semilogx()
    
def AddToPlot(times, lambdas, id=1):
    plt.figure(id)
    plt.step(times, lambdas, alpha=0.7)
    
def SavePlot(fout, id=1):
    plt.figure(id)
    plt.savefig(fout)
    
def PrintJAFSFile(jaf):
    print("#Migration_JAFS version 0.1")
    norm = sum(jaf)
    print("total\t", norm)
    #jaf = [v/norm for v in jaf]
    jfn = ["0100", "1100", "0001", "0101", "1101", "0011", "0111"]
    for v in zip(jaf, jfn):
        print(v[1], "\t", v[0])