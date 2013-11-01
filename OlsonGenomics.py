from OlsonUtilities import *

def chopSimple(fi, pathOut, size):
    '''Takes a file fi, chops it into as many segments of length size as 
    possible, and deposits the files containing the new segments in pathOut.
    Assumes only one genome per input file.'''
    
    f = open(fi,'r')
    ensureDir(pathOut+"/")
    nm = pathOut + "/" + fi.split("/")[-1][:-4]

    ctbp = 0
    ctsq = 0
    
    fout = open(nm+'_chop_{:03}kbp_{:02d}.fna'.format(int(size/1000),ctsq),'w')
    
    buf = f.readline()
    infostr = buf.rstrip()
    fout.write('{!s}; Seq # {!s}\n'.format(infostr,ctsq))
    buf = f.readline().rstrip()
    while buf:
        if ctbp <= size:
            fout.write(buf + "\n")
            ctbp += len(buf)
        else:
            fout.close()
            ctsq += 1
            ctbp = 0
            fout = open(nm+'_chop_{:03}kbp_{:02d}.fna'.\
                format(int(size/1000),ctsq),'w')
            fout.write('{!s}; Seq # {!s}\n'.format(infostr,ctsq))
            fout.write(buf + "\n")
            ctbp += len(buf)
        buf = f.readline().rstrip()
    fout.close()
    f.close()


def chopFolder(pathWork, sizeChop, numSeq):
    import os
    
    if len(sizeChop) != len(numSeq):
        numSeq = numSeq[0] * len(sizeChop)
    
    if pathWork[-1] != '/':
        pathWork = pathWork + '/'
    
    allList = os.listdir(pathWork)
    genomeList = []
    for file in allList:
        if file[-4:] == ".fna":
            genomeList.append(file)
    
    for file in genomeList:
        # make sequences to be matched
        for i in range(len(sizeChop)):
            sizestr = sizeString(sizeChop[i])
            ensureDir(pathWork+"Sequences_{!s}/".format(sizestr))
            #chopRandom(pathWork+file, pathWork+"Sequences_{!s}/".\
            #    format(sizestr), sizeChop[i], sizeChop[i]/10,\
            #        numSeq[i])
            chopRandomNoOverlap(pathWork+file, pathWork+"Sequences_{!s}/".\
                format(sizestr), sizeChop[i], sizeChop[i]/10,\
                    numSeq[i])
                

def analyzeChop(path, nm, out):
    import os
    import numpy
    dirList = os.listdir(path)
    listMatch = []
    nmlen = len(nm)
    for fname in dirList:
        if len(fname) > nmlen and fname[:nmlen] == nm and fname[-8:] == '_out.txt':
            listMatch.append(fname)

    klDivArr = numpy.zeros(len(listMatch), dtype = numpy.dtype(float))
    fracArr = numpy.zeros(len(listMatch), dtype = numpy.dtype(float))
    kmerDict = {}
    sz = 0

    for match in listMatch:
        v = int(match[-14:-12])-1
        f = open(path+"/"+match,'r')
        
        buf = f.readline().rstrip()
        sz = int(buf.split('\t')[1])

        buf = f.readline().rstrip()
        klDivArr[v] = float(buf.split('\t')[1])

        buf = f.readline().rstrip()
        fracArr[v] = float(buf.split('\t')[1])

        f.readline()
        buf = f.readline().rstrip()
        while len(buf) > 7:
            sp = buf.split()
            kmer = sp[0]
            ct = int(sp[1])
            if kmer in kmerDict.keys():
                kmerDict[kmer][v] = float(1000*ct)/float(sz)
            else:
                kmerDict[kmer] = numpy.zeros(len(listMatch),\
                    dtype = numpy.dtype(float))
                kmerDict[kmer][v] = float(1000*ct)/float(sz)
            buf = f.readline().rstrip()
        print v+1
        f.close()

    fout = open(path+"/"+out,'w')
    fout.write("Name: {0!s}\n\n".format(nm));
    fout.write("Fractions:\n\tMean = {:1.6f}\n\tStandard Deviation={:1.6f}\n\n".\
        format(numpy.mean(fracArr),numpy.std(fracArr)))
    fout.write("Kmers:\n")
    for kmer in sorted(kmerDict):
        arr = kmerDict[kmer]
        spread = numpy.amax(arr)-numpy.amin(arr)
        fout.write("\t{!s}: Mean = {:1.6f}, StdDev = {:1.6f}, Spread = {:1.6f}, Spread/Mean = {:1.6f}\n".format(kmer, numpy.mean(arr), numpy.std(arr), spread, spread/numpy.mean(arr)))
    fout.close()


def chopRandom(fi, pathOut, avg_size, margin, num_chop):
    import random, string
    _, seq = readSequence(fi)
    llim = len(seq) - avg_size - margin
    if llim <= 0:
        return
        
    sizestr = sizeString(avg_size)

    if pathOut[-1] != "/":
        pathOut = pathOut + "/"
    
    name = fi.split("/")[-1][:-4]
    for i in range(num_chop):
        f = open("{!s}{!s}_random_chopped_{!s}b_{!s}.fna".\
            format(pathOut,name,sizestr,str(i).zfill(3)),'w')
        st = random.randrange(llim)
        leng = random.randrange(avg_size - margin, avg_size + margin)
        subseq = seq[st:st+leng]
        f.write(">: {!s}, {!s}bp: {!s} - {!s}\n".\
            format(string.replace(name,"_"," "), sizestr, st, st+leng))
        #f.write(">: {!s}\n".format(name))
        ct = 0
        while 1:
            f.write(subseq[ct:ct+60] + "\n")
            ct += 60
            if (ct > len(subseq)):
                break
        f.close()


def chopRandomNoOverlap(fi, pathOut, avg_size, max_margin, num_chop):
    import random, string
    _, seq = readSequence(fi)
    interval = len(seq)/num_chop
    if interval < avg_size:
        return
    else:
        margin = min((interval-avg_size)/2, max_margin)
        
    sizestr = sizeString(avg_size)
        
    if pathOut[-1] != "/":
        pathOut = pathOut + "/"
    
    name = fi.split("/")[-1][:-4]
    for i in range(num_chop):
        f = open("{!s}{!s}_random_chopped_{!s}b_{!s}.fna".\
            format(pathOut,name,sizestr,str(i).zfill(3)),'w')
        st = random.randrange(margin)+(i*interval)
        leng = random.randrange(avg_size - margin, avg_size + margin)
        subseq = seq[st:st+leng]
        f.write(">: {!s}, ~{!s}bp: {!s} - {!s}\n".\
            format(string.replace(name,"_"," "), sizestr, st, st+leng))
        ct = 0
        while 1:
            f.write(subseq[ct:ct+60] + "\n")
            ct += 60
            if (ct > len(subseq)):
                break
        f.close()
        

def testCorrectPlus(raifi):
    import string
    f = open(raifi,'r')
    total = 0
    correct = 0
    reps = {}
    f.readline()
    buf = f.readline().rstrip()
    while buf:
        genus = buf.split(" | ")[-1].split(" ")[0]
        if genus not in reps:
            reps[genus] = []
        buf = f.readline().rstrip()
        if string.find(buf, genus) > -1:
            correct += 1
            sqno = buf.split("; Seq # ")[1]
            reps[genus].append(sqno)
        total += 1
        buf = f.readline().rstrip()
    return float(100.0 * correct)/float(total), reps


def testCorrectTopList(resFile, outFile):
    import string
    f = open(resFile, 'r')
    totalCount = 0
    innerCount = 0
    correctScore = [float('NaN')]
    correctRank = [-1]
    realGenus = ''
    line = []
    f.readline()
    buf = f.readline()
    while buf:
        if buf == "\n":
            totalCount += 1
        elif buf[0] == ">":
            realGenus = buf.split(":")[1].split("_")[0].lstrip()
            correctScore.append(float('NaN'))
            correctRank.append(-1)
            innerCount = 0
        else:
            line = buf.rstrip().split("\t")
            innerCount += 1
            if string.find(line[1], realGenus) > -1:
                correctScore[totalCount] = line[0]
                correctRank[totalCount] = innerCount
        buf = f.readline()
    f.close()

    oF = open(outFile, 'w')
    oF.write("Rank,Score\n")
    for l in zip(correctRank, correctScore):
        oF.write("{!s},{!s}\n".format(l[0],l[1]))
    oF.close()


def routeSeqs(family, pathIn, mainNo, mainPath, divertNo, otherPath):
    import os
    import shutil
    import string
    import random
    
    files = os.listdir(pathIn)
    fm = []
    for f in files:
        if string.find(f,family) > -1:
            fm.append(f)

    fms = set(fm)

    ensureDir(mainPath+"/")
    ensureDir(otherPath+"/")

    divertSet = set(random.sample(fms, divertNo))
    if mainNo + divertNo < len(fm):
        mainSet = random.sample(fms - divertSet, mainNo)
    else:
        mainSet = fms - divertSet

    for m in fm:
        if m in divertSet:
            shutil.move(pathIn + '/' + m, otherPath + '/' + m)
        elif m in mainSet:
            shutil.move(pathIn + '/' + m, mainPath + '/' + m)
        else:
            os.remove(pathIn + '/' + m)


def makeDistinct(llist, nums, reps, destPth, genePth, numSq, numDB, nKBP):
    import random
    for n in nums:
        for r in range(reps):
            subset = random.sample(llist, n)
            super = '{!s}{!s}sp_{!s}kbp_{!s}'.format(destPth, n, nKBP, r)
            for s in subset:
                fam = s.split('_')[0]
                mnpth = '{!s}/testseqs'.format(super)
                orpth = '{!s}/DB_input'.format(super)
                chopSimple('{!s}{!s}'.format(genePth, s), super, nKBP*1000)
                routeSeqs(fam, super, numSq, mnpth, numDB, orpth)


def testCorrectSubfolders(path):
    import os
    if path[-1] != '/':
        path = path + '/'
    li = os.listdir(path)
    corr = []
    for fl in li:
        res = path + fl + '/res'
        try:
            tc, rs = testCorrect(res)
            corr.append(tc)
        except IOError:
            print "{!s}: Whoops".format(fl)
    return corr


def randomSeedClustering(fi, pr, reps, pathWork, pathRai):
    import os
    import numpy
    import networkx as nx
    import community
    import matplotlib.pyplot as plt
    
    contigs = []
    input = open(fi,'r')
    buf = input.readline()
    while buf:
        contigs.append(buf.split('\t')[0])
        buf = input.readline()
    input.close()

    if pathWork[-1] != '/':
        pathWork = pathWork + '/'
    if pathRai[-1] != '/':
        pathRai = pathRai + '/'

    mat = numpy.zeros((len(contigs),len(contigs)), dtype=numpy.float)

    for i in range(reps):
        os.system("perl makeRandom.pl {!s} {!s} {!s}dbIn.fna {!s}seqs.fna".format(pr, fi, pathWork, pathWork))
        os.system("{!s}raiphy -e .fna -m 2 -i {!s}dbIn.fna -d {!s}db".format(pathRai, pathWork, pathWork))
        os.system("{!s}raiphy -e .fna -m 0 -i {!s}seqs.fna -d {!s}db -o {!s}res".format(pathRai, pathWork, pathWork, pathWork))
        os.system("rm {!s}dbIn.fna".format(pathWork))
        os.system("rm {!s}seqs.fna".format(pathWork))
        os.system("rm {!s}db".format(pathWork))
        
        res = open("{!s}res".format(pathWork))
        lines = res.readlines()
        for sq, db in zip(lines[1::2], lines[2::2]):
            a = sq.rstrip()[1:]
            b = db.rstrip()
            mat[contigs.index(a)][contigs.index(b)] += 1.0
        res.close()
    os.system("rm {!s}res".format(pathWork))

    mat = mat/reps # Normalizing

    adj = open("{!s}adjacency_{!s}.csv".format(pathWork, reps), 'w')
    for j in range(len(contigs)):
        adj.write(",".join(str(n) for n in mat[j])+"\n")
    adj.close()

    D = nx.Graph() # the community method requires an undirected graph
    D.add_nodes_from(contigs)
    for j in range(len(contigs)):
        x = contigs[j]
        for k in range(j):
            y = contigs[k]
            D.add_edge(x,y,weight=max(mat[j][k],mat[k][j]))

    gra = open("{!s}graph_facts_{!s}.txt".format(pathWork, reps),'w')
    gra.write("Communities: \n")
    # http://perso.crans.org/aynaud/communities/
    partition = community.best_partition(D)
    gra.write("Modularity: {!s}\n".format(community.modularity(partition, D)))
    size = float(len(set(partition.values())))
    pos = nx.spring_layout(D)
    count = 0.0
    for com in set(partition.values()):
        count = count + 1.0
        list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
        nx.draw_networkx_nodes(D, pos, list_nodes, node_size = 20,
                       node_color = [count/size]*len(list_nodes), with_labels=False,cmap=plt.cm.Dark2,vmin=0.0, vmax=1.0)
        gra.write("{!s}: {!s}\n".format(com,sorted(list_nodes)))
    nx.draw_networkx_edges(D, pos)
    plt.savefig("{!s}graph_{!s}.png".format(pathWork, reps))
    gra.close()


def bootStrap(fi, threshold, pathWork, pathRai):
    import os
    #import numpy
    #import scipy

    #contigs = []
    #input = open(fi,'r')
    #buf = input.readline()
    #while buf:
    #    contigs.append(buf.split('\t')[0])
    #    buf = input.readline()
    #input.close()

    #mat = numpy.zeros((len(contigs),len(contigs)), dtype=numpy.float)

    if pathWork[-1] != '/':
        pathWork = pathWork + '/'
    if pathRai[-1] != '/':
        pathRai = pathRai + '/'
    
    os.system("perl separateBySize.pl {!s} {!s} {!s}dbIn.fna {!s}seqs.fna {!s}nums.txt".format(threshold, fi, pathWork, pathWork, pathWork))
    nmf = open("{!s}nums.txt".format(pathWork),'r')
    [dbNum, seqNum] = nmf.readline().rstrip().split(",")
    os.system("{!s}raiphy -e .fna -m 2 -i {!s}dbIn.fna -d {!s}db".format(pathRai, pathWork, pathWork))
    os.system("{!s}raiphy -e .fna -m 0 -i {!s}seqs.fna -d {!s}db -o {!s}resSeeds".format(pathRai, pathWork, pathWork, pathWork))
    os.system("rm {!s}dbIn.fna".format(pathWork))
    os.system("rm {!s}seqs.fna".format(pathWork))
    os.system("rm {!s}db".format(pathWork))
    os.system("rm {!s}nums.txt".format(pathWork))

    seqs = []
    matches = []
    res = open("{!s}resSeeds".format(pathWork),'r')
    lines = res.readlines()
    for sq, db in zip(lines[1::2], lines[2::2]):
        seqs.append(int(sq.rstrip()[2:]))
        matches.append(int(db.rstrip()[1:]))
    seeds = sorted(set(matches))

    dict1 = {}
    for j in range(len(seeds)):
        dict1[seeds[j]] = [seeds[j]]
    for i in range(len(seqs)):
        s = matches[i]
        dict1[s].append(seqs[i])
    list1 = dict1.values()

    fout = open("{!s}Clusters_{!s}.txt".format(pathWork, threshold),'w')
    fout.write(str(list1))
    fout.close()


def adjMatrix2Dendrograms(fi, namefi, pathOut, num):
    import numpy
    import scipy.cluster.hierarchy as sch
    import pickle
    import pylab
    
    if pathOut[-1] != "/":
        pathOut = pathOut + "/"
    
    pathOut = pathOut + str(num)
    
    f = open(fi,'r')
    lines = f.readlines()
    mat = numpy.zeros((len(lines),len(lines)), dtype=numpy.float)

    for i in range(len(lines)):
        row = lines[i].rstrip().split(',')
        for j in range(len(lines)):
            mat[i][j] = float(row[j])
    f.close()

    D = numpy.zeros((len(lines),len(lines)), dtype=numpy.float)

    for i in range(len(lines)):
        for j in range(i):
            D[i][j] = mat[i][j] + mat[j][i] # make D a symmetric distance matrix
            D[j][i] = mat[i][j] + mat[j][i]

    nf = open(namefi,'r')
    nlines = nf.readlines()
    names = []
    for i in nlines:
        names.append(i.rstrip())
    nf.close()

# http://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python/3011894#3011894

    mx = numpy.max(D)

    fig = pylab.figure(figsize=(8,6))
    ax1 = fig.add_axes([0.09,0.1,0.2,0.8])
    Y = sch.linkage(D, method='centroid')
    Z1 = sch.dendrogram(Y, orientation='right')
    ax1.set_xticks([])
    ax1.set_yticks([])

    axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
    idx1 = Z1['leaves']
    mat = mat[idx1,:]
    idx1.reverse()
    mat = mat[:,idx1]
    im = axmatrix.matshow(mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
    pylab.colorbar(im, cax=axcolor)
    fig.savefig("{!s}_dendrogram_matrix.png".format(pathOut))

    order = open("{!s}_dendro_order.p".format(pathOut),'wb')
    pickle.dump(Z1['leaves'],order)
    order.close()

    clustersOut = open("{!s}_dendro_clusters.txt".format(pathOut),'w')

    for x in range(min(10,int(mx*20)),0,-1):
        clusters = sch.fcluster(Y, float(x)/10.0, criterion='distance')
        clustersOut.write("Distance threshold of {!s}\n".format(float(x)/10.0))
        mx = numpy.max(clusters)
        clustersOut.write("Number of clusters: {!s}\n".format(mx))
        clDict = {}
        for i in range(len(names)):
            if clusters[i] in clDict:
                clDict[clusters[i]].append(int(names[i]))
            else:
                clDict[clusters[i]] = [int(names[i])]
        for i in range(1,mx+1):
            clustersOut.write("Cluster #{!s}: {!s}\n".format(i,clDict[i]))
        clustersOut.write("\n")
    clustersOut.close()


def rai2KMeans(fi, clRange, out):
    import scipy.cluster.vq as scv
    import numpy
    
    names = []
    obs = []
    
    f = open(fi,'r')
    f.readline()
    buf = f.readline().rstrip()
    while buf:
        l = buf.split(":")
        names.append(l[0][1:])
        o = [0.0]*(len(l)-1)
        for i in range(len(l)-1):
            o[i] = float(l[i+1])
        obs.append(o)
                
        buf = f.readline().rstrip()

    D = numpy.array(obs)

    W = scv.whiten(D)
            
    clustersOut = open(out,'w')
    
    for k in clRange:
        centroids,_ = scv.kmeans(W, k)
        clusters,_ = scv.vq(W,centroids)
        clDict = {}
            
        for i in range(len(names)):
            if clusters[i] in clDict:
                clDict[clusters[i]].append(int(names[i]))
            else:
                clDict[clusters[i]] = [int(names[i])]
                    
        clustersOut.write("k: {!s}\n[".format(k))
        clustersOut.write(",".join(str(v) for v in clDict.values()))
        clustersOut.write("]\n\n")
    clustersOut.close()


def distMatrixInputOrder(matfile, orderfile, outfile):
    import pickle
    import numpy
    import pylab
    
    mf = open(matfile,'r')
    lines = mf.readlines()
    mat = numpy.zeros((len(lines),len(lines)), dtype=numpy.float)
        
    for i in range(len(lines)):
        row = lines[i].rstrip().split(',')
        for j in range(len(lines)):
            mat[i][j] = float(row[j])
    mf.close()

    order = open(orderfile,'rb')
    idx = pickle.load(order)
    order.close()

    fig = pylab.figure(figsize=(6,6))
    axmatrix = fig.add_axes([0.05,0.1,0.8,0.8])
    mat = mat[idx,:]
    idx.reverse()
    mat = mat[:,idx]
    im = axmatrix.matshow(mat, aspect='auto', origin='lower', cmap=pylab.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    axcolor = fig.add_axes([0.88,0.1,0.02,0.8])
    pylab.colorbar(im, cax=axcolor)
    fig.savefig(outfile)


def csv2DistanceDistribution(file, out, intervals = [], include = 1):
    # http://stackoverflow.com/a/5328669
    import numpy as np
    
    f = open(file,'r')
    lines = f.readlines()
    mylist = []
    nl = len(lines)
    if intervals == []:
        for i in range(nl):
            row = lines[i].rstrip().split(',')
            for j in range(nl):
                if float(row[j]) > 0:
                    mylist.append(float(row[j]))
    else:
        if include:
            for n in range(len(intervals)):
                for i in intervals[n]:
                    row = lines[i].rstrip().split(',')
                    for j in set(intervals[n]) - set([i]):
                        mylist.append(float(row[j]))
        else:
            for n in range(len(intervals)):
                for i in intervals[n]:
                    row = lines[i].rstrip().split(',')
                    okset = set(range(nl)) - set(intervals[n])
                    for j in okset:
                        mylist.append(float(row[j]))
    f.close()

    mat = np.array(mylist)
    np.savetxt(out, mat, fmt='%.5e', delimiter=',')


def getDistances(fi, pctrep):
    import numpy
    import random
        
    names = []
    obs = []
    
    f = open(fi,'r')
    f.readline()
    buf = f.readline().rstrip()
    while buf:
        l = buf.split(":")
        names.append(l[0][1:])
        o = [0.0]*(len(l)-1)
        for i in range(len(l)-1):
            o[i] = float(l[i+1])
        obs.append(o)
        
        buf = f.readline().rstrip()
    
    D = numpy.array(obs)

    numrep = int(pctrep * len(D)/100.0)
            
    reps = random.sample(xrange(len(D)),numrep)

    dists = []
    for i in range(len(reps)):
        for j in range(i):
            dists.append(numpy.linalg.norm(D[i]-D[j]))

    return dists


def distanceHistogramFolder(folder, pctrep, max):
    import matplotlib.pyplot as plt
    import os
    
    if folder[-1] != '/':
        folder = folder + '/'
    li = os.listdir(folder)
    dist = []
    for fl in li:
        if fl[-3:] == ".db" and len(dist) < max:
            dist.append(getDistances(folder + fl, pctrep))

    dist = [item for sublist in dist for item in sublist]

    # http://matplotlib.org/examples/api/histogram_demo.html

    plt.figure()
    n, bins, patches = plt.hist(dist, 50, facecolor='green')

    plt.savefig("{!s}Histogram.png".format(folder))


def speciesDistances(pathWork, pathRai, sizeChop, numSeq, type='euclidean'):
    import os
    
    if len(sizeChop) != len(numSeq):
        numSeq = numSeq[0] * len(sizeChop)
        
    if pathWork[-1] != '/':
        pathWork = pathWork + '/'
    if pathRai[-1] != '/':
        pathRai = pathRai + '/'
        
    allList = os.listdir(pathWork)
    genomeList = []
    for file in allList:
        if file[-4:] == ".fna":
            genomeList.append(file)

    # ns = len(genomeList)
        
    for file in genomeList:
        # make sequences to be matched
        ensureDir(pathWork+"Sequences/")
        for i in range(len(sizeChop)):
            chopRandom(pathWork+file, pathWork+"Sequences/", sizeChop[i], sizeChop[i]/10, numSeq[i])
        
    # Make RAI databases
    os.system("{!s}raiphy -e .fna -m 2 -I {!s}Sequences/ -d {!s}seqs".format(pathRai,pathWork,pathWork))
    
    computeDistances(pathWork+"seqs", pathWork+"dist_{!s}bp_{!s}reps_{!s}.csv".format(sizeChop,numSeq,type), method=type)

    os.system("rm -r {!s}Sequences/".format(pathWork))
    os.system("rm -r {!s}seqs".format(pathWork))

    print "Files removed"


def testSelfRecruitment(pathWork, pathRai, numSpecies, sizeChop, numDB, numSeq):
    import os
    import random
    import numpy as np
    import scipy.cluster.vq as scv
    import scipy.cluster.hierarchy as sch
    #import networkx as nx
    #import community as cm
    #import matplotlib.pyplot as plt

    if pathWork[-1] != '/':
        pathWork = pathWork + '/'
    if pathRai[-1] != '/':
        pathRai = pathRai + '/'
            
    allList = os.listdir(pathWork)
    genomeList = []
    for file in allList:
        if file[-4:] == ".fna":
            genomeList.append(file)

    subset = random.sample(genomeList, numSpecies)

    for file in subset:
        # make sequences to be matched
        ensureDir(pathWork+"Sequences/")
        chopRandom(pathWork+file, pathWork+"Sequences/", sizeChop, sizeChop/10, numSeq)
        # make sequences to be matched to
        ensureDir(pathWork+"DataBase/")
        chopRandom(pathWork+file, pathWork+"DataBase/", sizeChop, sizeChop/10, numDB)

    # Make RAI databases
    os.system("{!s}raiphy -e .fna -m 2 -I {!s}Sequences/ -d {!s}seqs".format(pathRai,pathWork,pathWork))
    os.system("{!s}raiphy -e .fna -m 2 -I {!s}DataBase/ -d {!s}db".format(pathRai,pathWork,pathWork))

    # Data sets for further evaluation
    namesSq, RaiSq = rai2Numpy(pathWork+"seqs")
    namesDb, RaiDb = rai2Numpy(pathWork+"db")
            
    # namesAll = namesSq + namesDb
    RaiAll = np.concatenate([RaiSq,RaiDb])

    csvout = open("{!s}seqs.csv".format(pathIn),'w')
    for r in RaiAll:
        csvout.write(",".join(str(x) for x in r)+"\n")
            
    # Run RAIphy
    os.system("{!s}raiphy -e .fna -m 0 -I {!s}Sequences/ -d {!s}db -o {!s}output".format(pathRai,pathWork,pathWork,pathWork))
    # Evaluate RAIphy results
    raiDict = {}
    for k in namesDb:
        raiDict[k] = [k]
    raiRes = open("{!s}output".format(pathWork),'r')
    raiRes.readline()
    buf = raiRes.readline().rstrip()
    while buf:
        key = buf[3:]
        buf = raiRes.readline().rstrip()
        sq = buf[1:]
        raiDict[key].append(sq)
        buf = raiRes.readline().rstrip()
    raiList = []
    for k in raiDict:
        raiList.append(raiDict[k])

    print "RAIphy cluster list: {!s}".format(raiList)

    os.system("rm -r {!s}Sequences/".format(pathWork))
    os.system("rm -r {!s}seqs".format(pathWork))
    os.system("rm -r {!s}DataBase/".format(pathWork))
    os.system("rm -r {!s}db".format(pathWork))
    os.system("rm -r {!s}output".format(pathWork))

    print "Files removed"

    # K-means
    whitened = newWhiten(RaiAll)
    centroidsNoSeeds, _ = scv.kmeans(whitened, numSpecies)
    idsNoSeeds, _ = scv.vq(whitened, centroidsNoSeeds)
    # Want to implement k-means with initial seeds, but holy moly things go wrong in spectacular fashion.
    # TODO: k-means with initial seeds.  Might have to code my own version.
  
    # Compute distance matrix and graph for other clusterings
    D = np.zeros((len(RaiAll),len(RaiAll)),dtype=np.float)
    for i in range(len(RaiAll)):
        for j in range(i):
            dist = np.linalg.norm(RaiAll[i]-RaiAll[j])
            D[i][j] = dist
            D[j][i] = dist
    #max = numpy.max(D)
    #G = nx.Graph()
    #for i in range(len(D)):
    #    for j in range(i):
    #        G.add_edge(namesAll[i], namesAll[j], weight = max-D[i][j])

    # hierarchy-based clusterings
    WeightedLink = sch.weighted(D)
    # mylen = len(namesDb)
    clustWeightedLink = sch.fcluster(WeightedLink,numSpecies,criterion='maxclust')
    hist, bins = np.histogram(clustWeightedLink, bins=numSpecies)
    print hist


def computeDistances(raiFile, outFile, intervals=[], method='euclidean', include=0):
    import numpy
    if method not in {'euclidean', 'negDist2', 'angle', 'similarityIndex','raiScore'}:
        print "Method {!s} not recognized.  Please try again."
    else:
        dispatch = {
                    'euclidean': numpy.linalg.norm,
                    'negDist2': negDist2,
                    'angle': spectralContrastAngle,
                    'similarityIndex': similarityIndex,
                    }
        _, arr = rai2Numpy(raiFile)
        dist = numpy.zeros([len(arr), len(arr)])
        if intervals == []:
            for i in range(len(arr)):
                for j in range(i):
                    d = dispatch[method](arr[i], arr[j])
                    dist[i][j] = d
                    dist[j][i] = d
        else:
            if include == 0:
                for I in intervals:
                    for i in I:
                        okset = set(range(len(arr))) - set(I)
                        for j in okset:
                            d = dispatch[method](arr[i],arr[j])
                            dist[i][j] = d
                            dist[j][i] = d
            else:
                for I in intervals:
                    for i in I:
                        for j in I:
                            d = dispatch[method](arr[i],arr[j])
                            dist[i][j] = d
                            dist[j][i] = d

        oF = open(outFile,'w')
        for a in dist:
            oF.write(",".join(str(x) for x in a)+"\n")
        oF.close()


def scoreDistribution(csvFile, outFile, intervals=[], include=0):
    import numpy
    scoreMat = numpy.genfromtxt(csvFile, dtype=numpy.float32, delimiter = ",")
    dist = []
    if include == 0:
        for I in intervals:
            for i in I:
                okset = set(range(len(scoreMat))) - set(I)
                for j in okset:
                    dist.append(scoreMat[i][j])
    else:
        for I in intervals:
            for i in I:
                for j in I:
                    if j != i:
                        dist.append(scoreMat[i][j])

    oF = open(outFile,'w')
    oF.write("\n".join(str(a) for a in dist))
    oF.close()


def correctDistributionOneAnswer(raiScoreFile, outFile, correctList):
    import numpy
    rsf = numpy.genfromtxt(raiScoreFile, dtype=numpy.float32, delimiter=",")
    nm = len(rsf)
    oF = open(outFile, 'w')
    oF.write("Correct:Max,Correct:Min,Rank(Correct)\n")
    for i in range(nm):
        ans = rsf[i][correctList[i]]
        sorted = numpy.sort(-rsf[i])
        sorted = -sorted
        oF.write("{!s},".format(ans/sorted[-1]))
        oF.write("{!s},".format(ans/sorted[0]))
        oF.write("{!s}\n".format(numpy.nonzero(sorted==ans)[0][0]+1))
    oF.close()


def JGI(fi, k):
    import matplotlib.pyplot as plt
    # initialize
    myFreqs = {}
    count = 0
    allFreqs = []
    allCounts = []
    allGCs = []
    
    # read file, get kmer frequencies and GC counts
    f = open(fi, 'r')
    buf = f.readline().rstrip()
    while buf:
        seq = ''
        buf = f.readline()
        while not buf.startswith('>') and buf:
            seq = seq + buf.rstrip()
            buf = f.readline()
        subFreq, subCount, subGC = kmerCountSequence(seq, k)
        count += subCount
        for km in subFreq:
            if km in myFreqs:
                myFreqs[km] += float(subFreq[km])
            else:
                myFreqs[km] = float(subFreq[km])
        allFreqs.append(subFreq)
        allCounts.append(subCount)
        allGCs.append(float(subGC)/float(subCount))
    for km in myFreqs:
        myFreqs[km] = myFreqs[km]/float(count)

    # compute y-axis values
    ratios = [0.0]*len(allCounts)
    for i in range(len(allCounts)):
        sum = 0.0
        for km in allFreqs[i]:
            sum += (float(allFreqs[i][km])/float(allCounts[i]))/myFreqs[km]
        ratios[i] = sum/float(4**k)

    # plot
    plt.plot(allGCs,ratios,'b.')
    plt.show()
