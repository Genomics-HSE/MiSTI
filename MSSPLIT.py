#!/Users/shchur/anaconda3/bin/python3
import sys

def WriteList(l, f, sep = ""):
    f.write(sep.join(l))
    f.write("\n")
    return
    for i in range(len(l)):
        f.write(l[i])
        if i < len(l) - 1:
            f.write(" ")
        else:
            f.write("\n")

if len(sys.argv) < 2:
    print("Need input file")
    exit(0)

fn = sys.argv[1]

jaf = [0 for i in range(7)]
fnw1 = "ms2g1.ms"
fnw2 = "ms2g2.ms"
fw1 = open(fnw1, 'w')
fw2 = open(fnw2, 'w')
#f = open(fn, 'r')
with open(fn) as f:
    for _ in range(2):        
        line = next(f)
        fw1.write(line)
        fw2.write(line)
    for line in f:        
        for _ in range(2):
            fw1.write(line)
            fw2.write(line)
            line = next(f)
        line = next(f)
        l = []
        line = line.split()
        l.append(line)
        line = next(f)
        l.append(line)
        line = next(f)
        l.append(line)
        line = next(f)
        l.append(line)
        line = next(f)
        l.append(line)
        f1 = [["positions:"],[],[]]
        f2 = [["positions:"],[],[]]
        for i in range( len(l[1])-1 ):
            if l[1][i] != l[2][i]:
                f1[0].append(l[0][i+1])
                f1[1].append(l[1][i])
                f1[2].append(l[2][i])
            if l[3][i] != l[4][i]:
                f2[0].append(l[0][i+1])
                f2[1].append(l[3][i])
                f2[2].append(l[4][i])
        fw1.write("segsites: " + str(len( f1[1] )) + "\n" )
        fw2.write("segsites: " + str(len( f2[1] )) + "\n"  )
        WriteList(f1[0], fw1, " ")
        WriteList(f2[0], fw2, " ")
        
        WriteList(f1[1], fw1, "")
        WriteList(f2[1], fw2, "")
        WriteList(f1[2], fw1, "")
        WriteList(f2[2], fw2, "")
fw1.close()
fw2.close()