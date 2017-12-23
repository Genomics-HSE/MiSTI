#!/Users/shchur/anaconda3/bin/python3
import sys

if len(sys.argv) < 2:
    print("Need input file")
    exit(0)

fn = sys.argv[1]

jaf = [0 for i in range(4)]
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
            s = int(l[0][i]) + int(l[1][i]) + int(l[2][i]) + int(l[3][i])
            #if s == 0:
            #    s = 4
            jaf[s-1] += 1
        i += 8
t = jaf[0] + jaf[1] + jaf[2]
print(jaf[0]/t)
print(jaf[1]/t)
print(jaf[2]/t)