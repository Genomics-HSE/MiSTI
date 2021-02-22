#!/usr/bin/env python3

#    Copyright (c) 2021 Vladimir Shchur (vlshchur@gmail.com)
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
import numpy as np
from scipy.optimize import least_squares
from math import exp

class PSMC:
    def __init__(self, inputFile, RD = -1, **kwargs):
        self.times = []
        self.eps = []
        self.pattern = []
        if inputFile:
            self.ReadPSMCFile(inputFile, RD)
        else:
            print("Currently PSMC initizalisation is possible from file only")

    def ReadPSMCFile(self, fn, RD = -1):
        maxRD = -1
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
        self.RD = RD

        with open(fn) as f:
            for line in f:
                line = line.split()
                if line[0] == "MM":
                    if line[1].startswith( "pattern" ):
                        pat = line[1][0:-1].split(":")[1]
                        pat = pat.split("+")
                        for p in pat:
                            p = [int(v) for v in p.split("*") ]
                            if len(p) == 1:
                                self.pattern.append( p[0] )
                            else:
                                self.pattern += [ p[1] for _ in range( p[0] )]
                    elif line[1].startswith( "n_seqs" ):
                        self.length = int(line[2][0:-1].split(":")[1])

                if line[0] != "RD" or int( line[1] ) != self.RD:
                    continue
                while line[0] != "RS":
                    if line[0] == "TR":
                        self.theta = float(line[1])
                        self.rho = float(line[2])
                    line = next(f)
                    line = line.split()
                while line[0] != "PA":
                    if line[0] != "RS":
                        print("Unexpected line.")
                        sys.exit(0)
                    self.times.append( float(line[2]) )
                    self.eps.append( float(line[3]) )
                    line = next(f)
                    line = line.split()
                break

    def ChangeTheta(self, theta):
        self.times = [v*self.theta/theta for v in self.times]
        self.eps = [v*self.theta/theta for v in self.eps]
        self.rho = self.rho/self.theta*theta
        self.theta = theta

    def CollapsePattern(self):
        times = []
        counter = 0
        for el in self.pattern:
            times.append( self.times[counter] )
            counter += el
        return( times )

    def AverageCoalescentRate(self, t1, t2):
        ci = 0
        self.times.append(np.inf)
        if t1 > t2:
            sys.exit(1)
        while self.times[ci] <= t1:
            ci += 1
        ci -= 1
        avCoalRate = 0.0
        et = 0.0
        tt = 0.0
        while ci < len(self.times)-1 and self.times[ci] < t2:
            tu = min( t2, self.times[ci+1] )
            tl = max( t1, self.times[ci] )
            ru = tu/self.eps[ci]
            rl = tl/self.eps[ci]
            avCoalRate += (ru-rl)
            tt += tu-tl
            ci += 1
        l1 = tt/avCoalRate
        self.times.pop(-1)
        return( l1 )

    def FitCoalescentTime(self, t1, t2):
        ci = 0
        self.times.append(np.inf)
        if t1 > t2:
            sys.exit(1)
        while self.times[ci] <= t1:
            ci += 1
        ci -= 1
        pnc = 0.0
        et = 0.0
        tt = 0.0
        while ci < len(self.times)-1 and self.times[ci] < t2:
            tu = min( t2, self.times[ci+1] )
            tl = max( t1, self.times[ci] )
            ru = tu/self.eps[ci]
            rl = tl/self.eps[ci]
            if ru == np.inf:
                vu = 0.0
            else:
                vu = (ru+1.0)*exp(rl-ru)
            et += exp(pnc)*((rl+1.0)-vu)*self.eps[ci]
            pnc -= (ru-rl)
            tt += tu-tl
            ci += 1
        et = et/( 1.0-exp(pnc) )
        l1 = least_squares(lambda l: (et - t1) - self.ExpectedCoalTime(l, tt), 1.0, bounds = (0.0, np.inf), ftol=4e-16, xtol=4e-16, gtol=4e-16)
        self.times.pop(-1)
        return( l1.x[0] )

    def ExpectedCoalTime(self, l, t):
        if t == np.inf:
            return(l)
        r = t/l
        v = (1.0-exp(-r)*(r+1.0))*l
        return( v/(1.0-exp(-r)) )

    def ReestimateCoalescentRates(self, times):
        t1 = times[0]
        ci = 0
        et = []
        for t1, t2 in zip(times[:-1],times[1:]):
            et.append( self.AverageCoalescentRate(t1, t2) )
        et.append(self.FitCoalescentTime(times[-1], np.inf)  )
        return(et)
