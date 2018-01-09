#!/Users/shchur/anaconda3/bin/python3
import sys

if len(sys.argv) < 2:
    print("Need input file")
    exit(0)

fn = sys.argv[1]

fnw1 = "data/ms2g1.ms"
fnw2 = "data/ms2g2.ms"
fw1 = open(fnw1, 'w')
fw2 = open(fnw2, 'w')

with open(fn) as f:
    for line in f:        
        for _ in range(2):
            fw1.write(line)
            fw2.write(line)
            line = next(f)
        chrLen = int( next(f) )
        lc = 0
        f1 = []
        f2 = []
        while 1:
            line = next(f)
            lc += 1
            if lc > int(chrLen):
                print("Too many segsites, expected at most " + str(chrLen))
                sys.exit(0)
            if line == "@end\n":
                break
            line = line.split("\t")
            if line[1][0] != line[1][1]:
                f1.append( line[0] )
            if line[1][2] != line[1][3]:
                f2.append( line[0] )
        fw1.write("@begin " + str( len(f1) ) + "\n" )
        fw1.write( chrLen + "\n" )
        for val in f1:
             fw1.write( val + "\t10" + "\n" )
        fw1.write("@end\n")
        
        fw2.write("@begin " + str( len(f2) ) + "\n" )
        fw2.write( chrLen + "\n" )
        for val in f2:
             fw2.write( val + "\t10" + "\n" )
        fw2.write("@end\n")
fw1.close()
fw2.close()