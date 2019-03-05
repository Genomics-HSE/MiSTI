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
import random

class MiPlot:#This is a class of static variables
    fig = None
    ax = None
    
class JAFS:#This is a class of static variables
    def __init__(self, jafs = [], ufLLHConst = None, fLLHConst = None, pop1 = None, pop2 = None):
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
        self.mi = None#migration rate
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
        if "mi" in kwargs:
            self.mi = kwargs["mi"]
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

def PrintErr(*args, sep="", endl="\n"):
    message = sep.join(args)
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
    if len(mu) == 0:
        llh = Migration.llh
    else:
        llh = Migration.JAFSLikelihood( mu )
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
        outData += "\t" + str(Migration.mi[i][0]) + "\t" + str(Migration.mi[i][1])
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
                    data.mi = [float(line[1]), float(line[2])]
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
        if data.mi is None:
            mu0_title = "-"
            mu1_title = "-"
        else:
            mu0_title = str(round(data.mi[0],1))
            mu1_title = str(round(data.mi[1],1))
        title = "llh = " + llh_title + ", migr (1->2) = " + mi1_title + ", migr (2->1) " + mi0_title + "\ninput file " + fmigr
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

def BootstrapJAFS(Jafs):
    genomeLen = 0
    for el in Jafs.jafs:
        if len(el) != 8:
            PrintErr("Cannot use provided SFS for bootstrap.")
            sys.exit(0)
        genomeLen += el[0]
    sfs = [0 for _ in range(8)]
    while sfs[0] < genomeLen:
        sfs_id = random.randint(0, len(Jafs.jafs)-1)
        for i in range(8):
            sfs[i] += Jafs.jafs[sfs_id][i]
    return(sfs)

def PrintJAFSFile(jaf, pop1 = False, pop2 = False):
    print("#MiSTI_JSFS version 1.0")
    if pop1:
        pop1 = pop1.strip("\n\r")
        print("#pop1", pop1, sep="\t")
    if pop2:
        pop2 = pop2.strip("\n\r")
        print("#pop2", pop2, sep="\t")
    jfn = ["total", "0100", "1100", "0001", "0101", "1101", "0011", "0111"]
    print("\t".join(jfn))
    if not isinstance(jaf, list):
        PrintErr("Unexpected SFS value: should be a list of a list of lists")
        sys.exit(0)
    if not isinstance(jaf[0], list):
        norm = sum(jaf)
        sfs_str = str(norm) + "\t" + "\t".join([str(v) for v in jaf])
        print(sfs_str)
    else:
        for sfs in jaf:
            if len(sfs) == 7:
                norm = sum(sfs)
                sfs_str = str(norm) + "\t" + "\t".join([str(v) for v in sfs])
                print(sfs_str)
            elif len(sfs) == 8:
                sfs_str = "\t".join([str(v) for v in sfs])
                print(sfs_str)
            else:
                print("Unexpected SFS entry.")
                sys.exit(0)
            

def ReadJAFS(fn, silent_mode=False):
    Jafs = JAFS()
    old_version = False
    with open(fn) as f:
        line = next(f, "EOF").rstrip("\n")
        if line.startswith("#MiSTI_JSFS") or line.startswith("#MiSTI_JAF") or line.startswith("#Migration_JAF"):
            pars = line.split(" ")
            if float(pars[2]) < 1:
                old_version = True
        else:
            PrintErr("Corrupted JSFS file header.")
            sys.exit(0)
    if old_version:
        return( ReadJAFS_old(fn, silent_mode) )
    with open(fn) as f:
        line = next(f, "EOF").rstrip("\n")
        if line.startswith("#MiSTI_JSFS"):
            pars = line.split(" ")
            if float(pars[2]) < 1:
                PrintErr("The file version is not supported anymore.")
                sys.exit(0)
        else:
            PrintErr("Corrupted JSFS file header.")
            sys.exit(0)
        while line[0] == "#":
            line = next(f, "EOF").rstrip("\n")
            if line[1:5] == "pop1":
                pars = line.split("\t")
                if len(pars) != 2:
                    PrintErr("Corrupted JSFS file header.")
                    sys.exit(0)
                Jafs.pop1 = pars[1]
                if not silent_mode:
                    print("pop1\t", pars[1])
            elif line[1:5] == "pop2":
                line = line.split("\t")
                if len(pars) != 2:
                    PrintErr("Corrupted JSFS file header.")
                    sys.exit(0)
                Jafs.pop2 = pars[1]
                if not silent_mode:
                    print("pop2\t", pars[1])
        if line.startswith("total"):
            line = next(f, "EOF").rstrip("\n")
        while line != "EOF":
            jsfs = line.split("\t")
            if len(jsfs) != 8:
                PrintErr("Unexpected line. Expected an entry for JSFS with eight TAB-separated columns.")
                sys.exit(0)
            Jafs.jafs.append([int(v) for v in jsfs])
            line = next(f, "EOF").rstrip("\n")
    return(Jafs)

def ReadJAFS_old(fn, silent_mode=False):
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
                if not silent_mode:
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
    Jafs.jafs.append(jafs)
    return(Jafs)    

def ReadMS(argument_string):
    PrintErr("WARNING: ReadMS() is not safe to use, the function has many assumptions on the ms command line")
    PrintErr("WARNING: assume that population sizes in the ancestral popualtion are equal for pop1 and pop2")
    #argument_string = "-n 2 3.0 -em 0.0 1 2 2.0 -em 0.05 2 1 3.0 -en 0.01 1 0.5 -en 0.02 2 0.05 -en 0.0375 1 0.5 -en 0.0375 2 0.5 -ej 1.25 2 1 -eM 1.25 0.0 -eN 1.25 1.0 -eN 2.0 5.0"
    args = argument_string.split(" ")
    pops = [[], []]
    migr = []
    puls = []
    splitT = 0
    popMerge = 0
    i = 0
    while i < len(args):
        if args[i] == "-n":
            pop = int(args[i+1])
            size = float(args[i+2])
            if pop != 1 and pop != 2:
                print("Population id should be 1 or 2.")
                print(args[i], args[i+1], args[i+2])
                sys.exit(0)
            pops[pop-1].append([0.0, size])
            i += 3
        elif args[i] == "-en":
            pop = int(args[i+2])
            time = float(args[i+1])
            size = float(args[i+3])
            if pop != 1 and pop != 2:
                print("Population id should be 1 or 2.")
                print(args[i], args[i+1], args[i+2], args[i+3])
                sys.exit(0)
            pops[pop-1].append([time, size])
            i += 4
        elif args[i] == "-eN":
            time = float(args[i+1])
            size = float(args[i+2])
            pops[0].append([time, size])
            pops[1].append([time, size])
            i += 3
        elif args[i] == "-em":
            time = float(args[i+1])
            direct = int(args[i+2])
            rate = float(args[i+4])
            migr.append([time, rate, direct])
            i += 5
        elif args[i] == "-es":#-es t i p
            time = float(args[i+1])
            pop = int(args[i+2])
            rate = 1 - float(args[i+3])
            puls.append([time, rate, pop])
            i += 4
        elif args[i] == "-ej":
            if int(args[i+2]) <= 2:
                splitT = float( args[i+1] )
                popMerge = int(args[i+2])
            i += 4
        else:
            i += 1
    if popMerge == 0:
        print("Populations should be merged. (-ej [time] 2 1)")
        sys.exit(0)
    if len(migr) > 2 or (len(migr) == 2 and migr[0][2] == migr[1][2]):
        print("Maximum two -em arguments are supported, for two arguments they should be in different directions")
        sys.exit(0)
    times = []
    for k in [0, 1]:
        for el in pops[k]:
            times.append(el[0])
    for el in migr:
        times.append(el[0])
    for el in puls:
        times.append(el[0])
    times.append(splitT)
    times = list(set(times))
    times.sort()
    i = 0
    while times[i] != splitT:
        i += 1 
    splitTind = i
    for k in [0,1]:
        pops[k].sort(key=lambda el: el[0])
        if pops[k][0][0] != 0.0:
            pops[k].insert(0, [0.0, 1.0])
    migr.sort(key=lambda el: el[0])
    puls.sort(key=lambda el: el[0])
    popSizes = [[0,0] for i in range(len(times))]
    for k in [0, 1]:
        ind = 0
        for i in range(len(times)):
            if times[i] == pops[k][ind][0]:
                popSizes[i][k] = pops[k][ind][1]
                curSize = pops[k][ind][1]
                ind += 1
            else:
                popSizes[i][k] = curSize
    mis = []
    ind = 0
    for i in range(splitTind):
        while ind < len(migr) and migr[ind][0] == times[i]:
            m = migr[ind]
            mis.append([m[2], i, splitTind, 2*m[1], 0])
            ind += 1
    pus = []
    ind = 0
    for i in range(splitTind):
        if ind == len(puls):
            break
        if times[i] == puls[ind][0]:
            p = puls[ind]
            pus.append([p[2], i, p[1], 0])
            ind += 1
    inputData = [None for _ in range(5)]
    inputData[0] = [2*(u-v) for u, v in zip(times[1:], times[:-1])]
    inputData[1] = [[1.0/u[0], 1.0/u[1]] for u in popSizes]
    inputData[2] = splitTind
    inputData[3] = mis#migration rates
    inputData[4] = pus#pulse migration rates
    return(inputData)

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
