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
import matplotlib.pyplot as plt
from math import log
from CorrectLambda import CorrectLambda
from MigrationInference import MigrationInference
import argparse

class MiPlot:#This is a class of static variables
    fig = None
    ax = None
    
class JAFS:#This is a class of static variables
    def __init__(self, jafs = None, ufLLHConst = None, fLLHConst = None, pop1 = None, pop2 = None):
        self.jafs = jafs
        self.ufLLHConst = ufLLHConst
        self.fLLHConst = fLLHConst
        self.pop1=pop1
        self.pop2=pop2
        
class MigData:
    def __init__(self, **kwargs):
        self.llh = None
        self.splitT = None
        self.migStart = None
        self.migEnd  = None
        self.times = None
        self.lambda1 = None
        self.lambda2 = None
        self.thrh = None
        self.mu = None#migration rate
        self.sampleDate = None
        if "llh" in kwargs:
            self.llh = kwargs["llh"]
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
        if "mu" in kwargs:
            self.mu = kwargs["mu"]
        if "sampleDate" in kwargs:
            self.sampleDate = kwargs["sampleDate"]

class Units:#This is a class of static variables
    mutRate = 1.25e-8
    binsize = 100
    N0 = 10000
    genTime = 1
    firstCall = True
    def __init__(self, **kwargs):
        if "mutRate" in kwargs:
            Units.muRate = kwargs["mutRate"]
        if "binsize" in kwargs:
            Units.binsize = kwargs["binsize"]
        if "N0" in kwargs:
            Units.N0 = kwargs["N0"]
        if "inpFile" in kwargs:
            Units.SetUnitsFromFile(kwargs["inpFile"])
        if Units.firstCall or len(kwargs) > 0:
            print("Units")
            print("mutation rate =", Units.mutRate, "\tbinsize =", Units.binsize, "\tN0 =", Units.N0, "\tgeneration time =", Units.genTime)
            Units.firstCall = False

    def SetUnitsFromFile(self, fn):
        try:
            with open(fn) as f:
                for line in f:
                    line = line.split("=")
                    if line[0] == "mutRate" and len(line) == 2:
                        try:#
                            Units.mutRate = float(line[1])
                        except:
                            print("Cannot read mutation rate entry from file, using default or previous values")
                            pass
                    elif line[0] == "binsize" and len(line) == 2:
                        try:
                            Units.binsize = float(line[1])
                        except:
                            print("Cannot read bin size entry from file, using default or previous values")
                            pass
                    elif line[0] == "N0" and len(line) == 2:
                        try:
                            Units.N0 = float(line[1])
                        except:
                            print("Cannot read N0 entry from file, using default or previous values")
                            pass
                    elif line[0] == "genTime" and len(line) == 2:
                        try:
                            Units.genTime = float(line[1])
                        except:
                            print("Cannot read generation time entry from file, using default values")
                            pass
        except:
            print("Units input file not found, using default values.")

def SetScaling():
    #mu = 1.1e-8
    print("\t\t\t\tPLEASE REMOVE ME!\n\t\t\tSINCERELY, SETSCALING()")
    u = Units()
    mutRate = 1.25e-8
    binsize = 100
    scaling = [mutRate, binsize]
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

def ReadPSMC(fn1, fn2, sampleDate = 0.0, RD = -1, doPlot = False):
    d1 = ReadPSMCFile(fn1, RD)
    d2 = ReadPSMCFile(fn2, RD)
#    if d1[2] != d2[2]:
#        print("Different RDs for input files 1 and 2.")
#        sys.exit(0)
    u = Units()
    scaleTime = d1[3]/(2.0*u.binsize*u.mutRate)*u.genTime
    scaleEPS = d1[3]/(2.0*u.binsize*u.mutRate)/2.0/u.N0
    
    d2[0] = [v*d2[3]/d1[3] for v in d2[0]]#rescale   time       by th1/th2
    d2[1] = [v*d2[3]/d1[3] for v in d2[1]]#rescale   epsize     by th2/th1 (compare with previous line!)
    
    sdResc = sampleDate/scaleTime
    if sdResc > 0:
        d2[0] = [v + sdResc for v in d2[0]]
        d2[0].insert(0, 0.0)
        d2[1].insert(0, 1.0)
    
    Tk = []
    Lk1 = []
    Lk2 = []
    Tk = sorted( d1[0] + d2[0][1:] )
    try:
        sampleDateDiscr = Tk.index(sdResc)
    except:
        print(Tk)
        print(sdResc)
        PrintErr("Unexpected error in ReadPSMC(). Get in touch with the author.")
        sys.exit(0)
    
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
    
    if doPlot:
        x = [v*scaleTime for v in Tk]
        y1 = [scaleEPS/v for v in Lk1]
        y2 = [scaleEPS/v for v in Lk2]
        AddToPlot(x, y1, "psmc1")
        AddToPlot(x[sampleDateDiscr:], y2[sampleDateDiscr:], "psmc2")
    L1tmp = [Lk1[i]]
    L2tmp = [Lk1[i]]
    Ttmp = [Tk[0]]
    Lk = [[u, v] for u, v in zip(Lk1, Lk2)]
    Tk = [ u - v for u, v in zip(Tk[1:], Tk[:-1])]
    return( [Tk, Lk, scaleTime, scaleEPS, d1[3], d1[4], sampleDateDiscr] )#time, coalescent rates, 2*N_0 (assuming default bin size = 100), effective population size/10000 rescale factor, theta and rho (from PSMC), sample date in discrite units

def OutputMigration(fout, mu, Migration):
    llh = Migration.JAFSLikelyhood( mu )
#    print( vars(Migration) )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times)+1)]   
    outData = "#MiSTI ver 0.3\n"
    outData += "LK\t" + str(llh) + "\n"#split time
    outData += "ST\t" + str(Migration.splitT) + "\n"#split time
    outData += "SD\t" + str(Migration.sampleDate) + "\n"#second sample date
    outData += "MS\t" + str(Migration.migStart) + "\n"#migration start
    outData += "ME\t" + str(Migration.migEnd) + "\n"#migration end
    outData += "MU\t" + str(mu[0]) + "\t" + str(mu[1]) + "\n"#migration
    outData += "TR\t" + str(Migration.thrh[0]) + "\t" + str(Migration.thrh[1]) + "\n"#migration
    outData += "SFS\t" + "\t".join(map(str, Migration.JAFS)) + "\n"#expected SFS
    dataJAFS = [v/sum(Migration.dataJAFS) for v in Migration.dataJAFS]
    outData += "DSF\t" + "\t".join(map(str, dataJAFS)) + "\n"#empirical SFS
    for i in range( len(times) ):
        outData += "RS\t" + str(times[i]) + "\t" + str(1.0/Migration.lc[i][0]) + "\t" + str(1.0/Migration.lc[i][1]) + "\n"
    
    if fout == "":
        print(outData)
    else:
        fw = open(fout, 'w')
        fw.write(outData)
        fw.close()

def OutputMigration2(fout, mu, Migration):
    if len(mu) == 0:
        llh = Migration.llh
    else:
        llh = Migration.JAFSLikelyhood( mu )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times)+1)]   
    outData = "#MiSTI2 ver 0.3\n"
    outData += "LK\t" + str(llh) + "\n"#split time
    outData += "ST\t" + str(Migration.splitT) + "\n"#split time
    outData += "SD\t" + str(Migration.sampleDate) + "\n"#second sample date
    outData += "TR\t" + str(Migration.thrh[0]) + "\t" + str(Migration.thrh[1]) + "\n"
    outData += "SFS\t" + "\t".join(map(str, Migration.JAFS)) + "\n"#expected SFS
    dataJAFS = [v/sum(Migration.dataJAFS) for v in Migration.dataJAFS]
    outData += "DSF\t" + "\t".join(map(str, dataJAFS)) + "\n"#empirical SFS
    for i in range( len(times) ):
        outData += "RS\t" + str(times[i]) + "\t" + str(1.0/Migration.lc[i][0]) + "\t" + str(1.0/Migration.lc[i][1])
        outData += "\t" + str(Migration.mu[i][0]) + "\t" + str(Migration.mu[i][1])
        if i < Migration.splitT:
            for val in Migration.Pr[i]:
                outData += "\t" + str(val[0]) + "\t" + str(val[1])
        outData += "\n"
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
    mu1 = []
    mu2 = []
    pr11 = [[],[]]
    pr22 = [[],[]]
    pr12 = [[],[]]
    sampleDate = 0
    data = MigData()
    with open(fmigr) as f:
        line = next(f).rstrip()
        line = line.split(" ")
        version = float(line[2])
        print("Format version: ", version)
        if version < 0.3:
            PrintErr("File version is not supported anymore.")
            sys.exit(0)
        if line[0] == "#MiSTI2":
            for line in f:
                line = line.split("\t")
                if line[0] == "LK":
                    data.llh = float(line[1])
                elif line[0] == "ST":
                    data.splitT = int(line[1])
                elif line[0] == "SD":
                    data.sampleDate = int(line[1])
                elif line[0] == "TR":
                    data.thrh = [float(line[1]), float(line[2])]
                elif line[0] == "SFS":
                    data.jaf = map(float, line[1:])
                elif line[0] == "RS":
                    times.append( float(line[1])*scaleTime )
                    lc1.append( 1.0/float(line[2])/scaleEPS )
                    lc2.append( 1.0/float(line[3])/scaleEPS )
                    mu1.append( float(line[4]) )
                    mu2.append( float(line[5]) )
                    if len(line) > 6:
                        pr11[0].append(float(line[6]))
                        pr11[1].append(float(line[7]))
                        pr22[0].append(float(line[8]))
                        pr22[1].append(float(line[9]))
                        pr12[0].append(float(line[10]))
                        pr12[1].append(float(line[11]))
                    else:
                        pr11[0].append(0)
                        pr11[1].append(0)
                        pr22[0].append(0)
                        pr22[1].append(0)
                        pr12[0].append(0)
                        pr12[1].append(0)
        else:
            for line in f:
                line = line.split("\t")
                if line[0] == "LK":
                    data.llh = float(line[1])
                elif line[0] == "ST":
                    data.splitT = int(line[1])
                elif line[0] == "SD":
                    data.sampleDate = int(line[1])
                elif line[0] == "MS":
                    data.migStart = int(line[1])
                elif line[0] == "ME":
                    data.migEnd = int(line[1])
                elif line[0] == "MU":
                    data.mu = [float(line[1]), float(line[2])]
                elif line[0] == "TR":
                    data.thrh = [float(line[1]), float(line[2])]
                elif line[0] == "SFS":
                    data.jaf = map(float, line[1:])
                elif line[0] == "RS":
                    times.append( float(line[1])*scaleTime )
                    lc1.append( 1.0/float(line[2])/scaleEPS )
                    lc2.append( 1.0/float(line[3])/scaleEPS )
    if doPlot:
        lc1 = [min(1.0/v,10) for v in lc1]
        lc2 = [min(1.0/v,10) for v in lc2]
#        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc1])
#        plt.step([v*scaleTime for v in times], [1.0/max(v,0.1)*scaleEPS for v in lc2])
        if data.llh == None:
            llh_title = "-"
        else:
            llh_title = str(round(data.llh,1))
        if data.mu is None:
            mu0_title = "-"
            mu1_title = "-"
        else:
            mu0_title = str(round(data.mu[0],1))
            mu1_title = str(round(data.mu[1],1))
        title = "llh = " + llh_title + ", migr (1->2) = " + mu1_title + ", migr (2->1) " + mu0_title + "\ninput file " + fmigr
        AddTitle(title)
        AddToPlot(times, lc1, "misti1")
        AddToPlot(times[sampleDate:], lc2[sampleDate:], "misti2")
        if len(pr11[0]) > 0:
            AddProb(pr11, pr22, pr12, times)
        splT=times[data.splitT]
        if data.migStart != None and data.migEnd != None:
            ms = times[data.migStart]
            me = times[data.migEnd]
            MiPlot.ax.axvspan(ms, me, color='k', alpha=0.05)
        MiPlot.ax.axvline(splT, color='k', alpha=0.1)
#    data = MigData(splitT = splitT, migStart = migStart, migEnd = migEnd, times = times, lambda1 = lc1, lambda2 = lc2, thrh = thrh, mu = mu, sampleDate = sampleDate, llh = llh)
    data.times = times
    data.lambda1 = lc1
    data.lambda2 = lc2
    return(data)

def ReadJAFS(fn, silent_mode=False):
    Jafs = JAFS()
    jafs = []
    with open(fn) as f:
        line = next(f).rstrip()
        while line[0] == "#":
            if line[1:10] == "MiSTI_JAF" or line[1:14] == "Migration_JAF":
                line = line.split(" ")
                if len(line) < 3:
                    PrintErr("Corrupted JAF file header.")
                    sys.exit(0)
                print("JAFS format version:", line[2])
            elif line[1:5] == "pop1":
                line = line.split(" ")
                if len(line) != 2:
                    PrintErr("Corrupted JAF file header.")
                    sys.exit(0)
                Jafs.pop1 = line[1]
                if not silent_mode:
                    print("pop1\t", line[1])
            elif line[1:5] == "pop2":
                line = line.split(" ")
                if len(line) != 2:
                    PrintErr("Corrupted JAF file header.")
                    sys.exit(0)
                Jafs.pop2 = line[1]
                if not silent_mode:
                    print("pop2\t", line[1])
            elif line[1:11] == "ufLLHConst":
                line = line.split(" ")
                if len(line) != 2:
                    PrintErr("Corrupted JAF file header.")
                    sys.exit(0)
                Jafs.ufLLHConst = float(line[1])
            elif line[1:10] == "fLLHConst":
                line = line.split(" ")
                if len(line) != 2:
                    PrintErr("Corrupted JAF file header.")
                    sys.exit(0)
                Jafs.fLLHConst = float(line[1])
            line = next(f).rstrip()

        line = line.split("\t")
        if len(line) != 2:
            PrintErr("Unexpected line. Expected an entry for JAFS with two TAB-separated columns.")
            sys.exit(0)
        jafs.append( int(line[1]) )
        for line in f:
            line = line.split("\t")
            if len(line) != 2:
                PrintErr("Unexpected line. Expected an entry for JAFS with two TAB-separated columns.")
                sys.exit(0)
            jafs.append( int(line[1]) )
    if len(jafs) != 8:
        print("Unexpected number of lines in the JAFS file.")
        sys.exit(0)
    Jafs.jafs = jafs
    return(Jafs)

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
    mutRate = scaling[0]
    binsize = scaling[1]
    scale = theta/(2.0*binsize*mutRate)
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
#    plt.figure(id)
    MiPlot.fig, (MiPlot.ax, MiPlot.pr11, MiPlot.pr22, MiPlot.pr12, MiPlot.nc) = plt.subplots(5, 1, gridspec_kw=dict(hspace=0.5, height_ratios=[3, 1, 1, 1, 1]))
    MiPlot.ax.semilogx()
    MiPlot.pr11.semilogx()
    MiPlot.pr22.semilogx()
    MiPlot.pr12.semilogx()
    MiPlot.nc.semilogx()
    
def AddTitle(title, id=1):
    MiPlot.ax.set_title(title)

def AddToPlot(times, lambdas, lbl = "", id=1):
    #plt.figure(id)
    MiPlot.ax.step(times+[2*times[-1]], [lambdas[0]]+lambdas, alpha=0.7, label=lbl)

def AddProb(pr11, pr22, pr12, times):
#    MiPlot.pr11 = plt.subplot(212)#, sharex = True)
    #, MiPlot.pr22, MiPlot.pr12, MiPlot.nc)
    nc = [None, None]
    nc[0] = [pr11[0][i]+pr22[0][i]+pr12[0][i] for i in range(len(pr11[0]))]
    nc[1] = [pr11[1][i]+pr22[1][i]+pr12[1][i] for i in range(len(pr11[1]))]
    for i in [0,1]:
        pr11[i] = [u/(v if v > 0 else 1) for u, v in zip(pr11[i], nc[i])]
        pr22[i] = [u/(v if v > 0 else 1) for u, v in zip(pr22[i], nc[i])]
        pr12[i] = [u/(v if v > 0 else 1) for u, v in zip(pr12[i], nc[i])]
    
    MiPlot.pr11.step(times+[2*times[-1]], [pr11[0][0]]+pr11[0], alpha=0.7, label="1")
    MiPlot.pr11.step(times+[2*times[-1]], [pr11[1][0]]+pr11[1], alpha=0.7, label="2")
    MiPlot.pr11.legend(loc="upper right", prop=dict(size=6))
    
    MiPlot.pr22.step(times+[2*times[-1]], [pr22[0][0]]+pr22[0], alpha=0.7, label="1")
    MiPlot.pr22.step(times+[2*times[-1]], [pr22[1][0]]+pr22[1], alpha=0.7, label="2")
    MiPlot.pr22.legend(loc="upper right", prop=dict(size=6))

    MiPlot.pr12.step(times+[2*times[-1]], [pr12[0][0]]+pr12[0], alpha=0.7, label="1")
    MiPlot.pr12.step(times+[2*times[-1]], [pr12[1][0]]+pr12[1], alpha=0.7, label="2")
    MiPlot.pr12.legend(loc="upper right", prop=dict(size=6))

    MiPlot.nc.step(times+[2*times[-1]], [nc[0][0]]+nc[0], alpha=0.7, label="1")
    MiPlot.nc.step(times+[2*times[-1]], [nc[1][0]]+nc[1], alpha=0.7, label="2")
    MiPlot.nc.legend(loc="upper right", prop=dict(size=6))

def SavePlot(fout, id=1):
    #plt.figure(id)
    MiPlot.ax.legend()
    MiPlot.fig.savefig(fout)
    
def PrintJAFSFile(jaf, pop1 = False, pop2 = False):
    print("#MiSTI_JAFS version 0.2")
    if pop1:
        pop1 = pop1.strip("\n\r")
        print("#pop1", pop1)
    if pop2:
        pop2 = pop2.strip("\n\r")
        print("#pop2", pop2)
    
    norm = sum(jaf)
    
    ufLLHConst = 0.0
    for j in range(1, norm+1):
        ufLLHConst += log(j)
    for i in range(7):
        for j in range(1, jaf[i]+1):
            ufLLHConst -= log(j)
    
    fLLHConst = 0.0
    for j in range(1, norm+1):
        fLLHConst += log(j)
    for j in range(1, jaf[0]+jaf[6]+1):
        fLLHConst -= log(j)
    for j in range(1, jaf[1]+jaf[5]+1):
        fLLHConst -= log(j)
    for j in range(1, jaf[2]+jaf[4]+1):
        fLLHConst -= log(j)
    for j in range(1, jaf[3]+1):
        fLLHConst -= log(j)
    
    print("#ufLLHConst", ufLLHConst)
    print("#fLLHConst", fLLHConst)
    
    print("total\t", norm)
    #jaf = [v/norm for v in jaf]
    jfn = ["0100", "1100", "0001", "0101", "1101", "0011", "0111"]
    for v in zip(jaf, jfn):
        print(v[1], "\t", v[0])




#Example of using ms parameters to run tests for perfectly known data
if 0:
    timesMS = [0, 0.0275, 0.0475, 0.175, 0.75, 3.75, 10]
    epsMS = [[13.0, 0.25], [0.5, 0.4], [0.5, 0.5], [3, 3], [2, 2], [3, 3], [6, 6]]
    inputData[0] = [2*(u-v) for u, v in zip(timesMS[1:], timesMS[:-1])]
    inputData[1] = [[1.0/u[0], 1.0/u[1]] for u in epsMS]

if 0:
    timesMS = [0, 0.0075, 0.0225, 0.125, 0.5, 1.875]
    epsMS = [[0.2, 0.2], [0.2, 0.45], [0.5, 0.5], [0.3, 0.3], [0.5, 0.5], [1.3, 1.3]]
    discr = 50
    timesTmp = [0]
    epsTmp = []
    for i in range(len(timesMS)-1):
        maxTime = timesTmp[-1]
        deltaT = timesMS[i+1] - timesMS[i]
        timesTmp += [maxTime + (j+1)*deltaT/discr for j in range(discr)]
        epsTmp += [epsMS[i] for j in range(discr)]
    epsTmp.append(epsMS[-1])
    timesMS = timesTmp
    epsMS = epsTmp
    inputData[0] = [2*(u-v) for u, v in zip(timesMS[1:], timesMS[:-1])]
    inputData[1] = [[1.0/u[0], 1.0/u[1]] for u in epsMS]
    inputData[2] = 20000
    inputData[3] = 1
