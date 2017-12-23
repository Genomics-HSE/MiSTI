#!/Users/shchur/anaconda3/bin/python3
import sys

if len(sys.argv) < 2:
    print("Need input file")
    exit(0)

fn = sys.argv[1]

jaf = [0 for i in range(7)]
#f = open(fn, 'r')
with open(fn) as f:
    for _ in range(2):        
        line = next(f)
    for line in f:
        for _ in range(4):
            line = next(f)
        l = []
        l.append(line)
        line = next(f)
        l.append(line)
        line = next(f)
        l.append(line)
        line = next(f)
        l.append(line)
        for i in range( len(l[0])-1 ):
            s0 = int(l[0][i]) + int(l[1][i])
            s1 = int(l[2][i]) + int(l[3][i])
            if s0 == 0:
                if s1 == 1:
                    jaf[2] += 1
                elif s1 == 2:
                    jaf[5] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
            if s0 == 1:
                if s1 == 0:
                    jaf[0] += 1
                elif s1 == 1:
                    jaf[3] += 1
                elif s1 == 2:
                    jaf[6] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
            if s0 == 2:
                if s1 == 0:
                    jaf[1] += 1
                elif s1 == 1:
                    jaf[4] += 1
                else:
                    self.PrintError("StateToJAF", "unexpected state " + self.PrintState(state) )
        i += 8
norm = sum(jaf)
print("-------------------",jaf[0]/norm,jaf[1]/norm,sep="\t\t")
print(jaf[2]/norm,jaf[3]/norm,jaf[4]/norm,sep="\t\t")
print(jaf[5]/norm,jaf[6]/norm,"-------------------",sep="\t\t")
print("singletons", (jaf[0]+jaf[2])/norm)
print("doubletons", (jaf[1]+jaf[3]+jaf[5])/norm)
print("tripletons", (jaf[4]+jaf[6])/norm)