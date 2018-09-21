lambdas, times, mu, splitT, theta = [], [], [], None, None

l = [1,1/0.4]
#l = [1]
lambdas = []
for la in l:
    lambdas.append([la,la])
#print(lambdas)
#times = [2*0.125]
#times = []
#lambdas = [[1,3],[0.1,0.5],[2,4],[1,0.5],[10,10]]
#times = [2,3,1,5]
#mu = [0.0, 0.0]#[0.00001,0.00001]
#mu = [0.00001,0.00001]
#splitT = 3

#times = [2*0.125]
lambdas = [[1.0, 2.5], [1.4285714285714286, 1.4285714285714286]]
mu = [100.0, 300.0]
mu = [0.0, 0.0]
times = [2.5]
splitT = 1

lambdas = [[1.0, 2.5], [1.4285714285714286, 1.6666666666666667], [2.0, 2.0]]
mu = [100.0, 30.0]
times = [0.25, 0.35]
splitT = 2
dataJAFS = [0 for i in range(7)]

inputData = [times, lambdas]







mu = [0.25, 0.25]
#migrR = inputData[1][0][0]/2
migrUnit = inputData[1][0][0]/2
mu = [migrUnit, migrUnit]
#mu = [0.0, 0.0]
'''theta = 0.000001
lambdas = [[1.0, 1.0], [2.5, 2.5]]
mu = [0.01, 0.003]
times = [0.25]
splitT = 0'''
if False:
    print("Input parameters:")
    print("\ttimes  = ", inputData[0], sep = "")
    print("\tlambda = ", inputData[1], sep = "")
    print("\tJAFS   = ", dataJAFS, sep = "")
    print("\tmigrat = ", mu, sep = "")
    print("\tsplitT = ", splitT, sep = "")
    print("\ttheta  = ", theta, sep = "")
    print("End of input.")
llh = []
if False:
    cl = CorrectLambda()
    mu1=[1/2,1/2]
    inter1 = [0, 1, 1, mu1[0], mu1[1]]
    inter2 = [0.02, 0.05, 0.05, mu1[0], mu1[1]]
    inter3 = [0.075, 0.5, 0.5, mu1[0], mu1[1]]
    inter4 = [2.5, 1, 1, mu1[0], mu1[1]]
    splitT = 2.5
    intervals = [inter1, inter2, inter3, inter4]
    coalRates = cl.CoalRates(intervals, splitT, 100)
    '''    print("./foreign/msHOT-lite/msHOT-lite 2 1000 -T -t 8196 -r 1355 3000000 -l ", end="")
    #print("./foreign/msHOT-lite/msHOT-lite 2 1000 -T -t 8196 -r 1355 3000000 -l -I 2 2 0 -em 0 1 2 1 -em 0 2 1 1 ", end="")
    for i in range(len(coalRates[0])):
        print(" -eN ", coalRates[1][i]/2.0, " ", 1/coalRates[0][i], end="", sep="")
    print(" -eN 1.25 1 > test_hmm/sim_coal_rate.ms")
    #print(" -eN 1.25 1 -eM 1.25 0 -ej 1.25 2 1 > data_std_migr_same1/sim_coal_rate.ms")'''
    coalRates[1] = coalRates[1] + [2.5, 5.0]
    coalRates[0] = coalRates[0] + [1, 1]
    plt.plot([v*2*10000 for v in coalRates[1]], [1.0/v for v in coalRates[0]])
    times = numpy.diff(coalRates[1])
    lambdas = [[v, v] for v in coalRates[0]]
    #plt.savefig("temp3.png")
    #sys.exit(0)
    cl.SetMu(mu1[0], mu1[1])
    print(times)
    print(len(lambdas))
    Migration = MigrationInference(times, lambdas, dataJAFS, mu1, 302, theta, correct = True, smooth = False)
    print("JAFSLikelyhood ideal=", exp(Migration.JAFSLikelyhood( mu1 )))
#    Migration.CorrectLambdas()
    plt.plot([v*2*10000 for v in coalRates[1][:len(Migration.lc)]], [1.0/v[0] for v in Migration.lc])
    if False:
        Migration = MigrationInference([0.02, 0.055, 2.5-0.075], [[1,1],[1/0.05,1/0.05],[1/0.5,1/0.5],[1,1]], dataJAFS, mu1, 3, theta, correct = False, smooth = False)
        print("JAFSLikelyhood ideal=", exp(Migration.JAFSLikelyhood( [1/4,1/4] )))
        print("JAFSLikelyhood ideal=", exp(Migration.JAFSLikelyhood( [1,1] )))
        print("JAFSLikelyhood ideal=", exp(Migration.JAFSLikelyhood( [1/2,1/2] )))
        print("JAFSLikelyhood ideal=", exp(Migration.JAFSLikelyhood( [4,4] )))
#    plt.savefig("temp3.png")
#    sys.exit(0)
#    print(exp( Migration.JAFSLikelyhood( mu1 ) ))
    
#    plt.savefig("temp3.png")
#for splitT in [20, 80, 90, 94, 95, 96, 97]:
for splitT in [98,99,100,101]:
    continue
    print("splitT = ", splitT)
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, correct = True, enableOutput = False, smooth = False)
#Migration.Solve()
    llh = Migration.MaxLLH(migrUnit)
    print("splitT = ", splitT, "\ttime = ", sum(inputData[0][0:splitT])*inputData[2], "\tllh = ", llh[0], "\tmu = ", [llh[1][0]/migrUnit,llh[1][1]/migrUnit])
#Migration.Solve()
#Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, correct = True, enableOutput = False, smooth = False)
#splitT = 97







for splitT in [80,85,90,95,102,103,104,105,106,107,108,109,110]:
    continue
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, correct = True, enableOutput = False, smooth = False)
    muSol = Migration.Solve()
    if muSol[1] > maxllh:
        maxmu = muSol[0]
        maxsplitT = splitT
    print(muSol)
    print("splitT = ", splitT, "\ttime = ", sum(inputData[0][0:splitT])*inputData[2], "\tmu = ", [muSol[0][0]/migrUnit,muSol[0][1]/migrUnit], "\tllh = ", muSol[1])

mu = [migrUnit, migrUnit]
print( Migration.JAFSLikelyhood( mu ) )

mu = [2*migrUnit, 2*migrUnit]
print( Migration.JAFSLikelyhood( mu ) )







for i in range(0, len(inputData[0]) - 3 ):
    if i < 80:
        llh.append(0)
        continue
#    if i != 95 and i != 96:
#         llh.append(0)
#         continue
    splitT = i
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splitT, theta, correct = True, enableOutput = False, smooth = False)
    llh.append( exp( Migration.JAFSLikelyhood( mu ) ) )
norm = sum(llh)
norm = 1.0
splT = [0,0,0]
for i in range( len(llh) ):
    if llh[i] > splT[1]:
        splT = [sum(inputData[0][0:i]), llh[i],i]
    print("SplitT\t", sum(inputData[0][0:i])*inputData[2], "(", i,")", "\tLLH\t", llh[i]/norm)
print("Maximum likelihood split time\t", splT[0]*inputData[2], "\tLLH\t", splT[1]/norm)
if doPlot:
    Migration = MigrationInference(inputData[0], inputData[1], dataJAFS, mu, splT[2], theta, correct = True, smooth = True)
    Migration.JAFSLikelyhood( mu )
    times = [sum(Migration.times[0:i]) for i in range(len(Migration.times))]
    plt.step([v*inputData[2] for v in times[0:]], [1.0/max(v[0],0)*inputData[3] for v in Migration.lc[0:-1]])
    plt.axvline(splT[0]*inputData[2], color='r')
    plt.savefig("temp3.png")

t2 = time.clock()
MigrationInference.Report()
print("Total time ", t2-t1)