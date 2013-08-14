def chopSimple(fi, pathout, size):
    '''Takes a file fi, chops it into as many segments of length size as possible, and deposits the files containing the new segments in pathout.
        Assumes only one genome per file.'''
    
    f = open(fi,'r')
    ensureDir(pathout+"/")
    nm = pathout + "/" + fi.split("/")[-1][:-4]

    ctbp = 0
    ctsq = 0
    
    fout = open(nm + '_chop_{:03}kbp_{:02d}.fna'.format(int(size/1000),ctsq),'w')
    
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
            fout = open(nm + '_chop_{:03}kbp_{:02d}.fna'.format(int(size/1000),ctsq),'w')
            fout.write('{!s}; Seq # {!s}\n'.format(infostr,ctsq))
            fout.write(buf + "\n")
            ctbp += len(buf)
        buf = f.readline().rstrip()
    fout.close()
    f.close()

def analyzeChop(path, nm, out):
    import os
    import string
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
                kmerDict[kmer] = numpy.zeros(len(listMatch), dtype = numpy.dtype(float))
                kmerDict[kmer][v] = float(1000*ct)/float(sz)
            buf = f.readline().rstrip()
        print v+1
        f.close()

    fout = open(path+"/"+out,'w')
    fout.write("Name: {0!s}\n\n".format(nm));
    fout.write("Fractions:\n\tMean = {:1.6f}\n\tStandard Deviation={:1.6f}\n\n".format(numpy.mean(fracArr),numpy.std(fracArr)))
    fout.write("Kmers:\n")
    for kmer in sorted(kmerDict):
        arr = kmerDict[kmer]
        spread = numpy.amax(arr)-numpy.amin(arr)
        fout.write("\t{!s}: Mean = {:1.6f}, StdDev = {:1.6f}, Spread = {:1.6f}, Spread/Mean = {:1.6f}\n".format(kmer, numpy.mean(arr), numpy.std(arr), spread, spread/numpy.mean(arr)))
    fout.close()

def chopRandom(fi, pathout, avg_size, interval, num_chop):
    import random
    nm, seq = readSequence(fi)
    llim = len(seq) - avg_size - interval
    if llim <= 0:
        return
    
    if avg_size >= 1000:
        sizestr = str(avg_size/1000) + "k"
    else:
        sizestr = str(avg_size)

    if pathout[-1] != "/":
        pathout = pathout + "/"
    
    out = fi[:-4]
    name = fi.split("/")[-1][:-4]
    f = open(pathout + name + "_random_chopped_" + sizestr + "b.fna",'w')
    for i in range(num_chop):
        st = random.randrange(llim)
        leng = random.randrange(avg_size - interval, avg_size + interval)
        subseq = seq[st:st+leng]
        #f.write(">: {!s}, {!s}bp: {!s} - {!s}\n".format(name, sizestr, st, st+leng))
        f.write(">: {!s}\n".format(name))
        ct = 0
        while 1:
            f.write(subseq[ct:ct+60] + "\n")
            ct += 60
            if (ct > len(subseq)):
                break
        f.write("\n")
    f.close()

def readSequence(fi):
    f = open(fi,'r')
    concat = ''
    buf = f.readline().rstrip()
    while buf:
        seq_name, seq = buf[2:], ''
        buf = f.readline()
        while not buf.startswith('>') and buf:
            seq = seq + buf.rstrip()
            buf = f.readline()
        if seq_name.find("complete genome") == -1:
            concat = concat + seq
        elif seq_name.find("complete genome") > -1:
            concat = seq
            break
    return seq_name, concat
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

def ensureDir(pth):
# source: http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    import os
    d = os.path.dirname(pth)
    if not os.path.exists(d):
        os.makedirs(d)

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

def randomSeedClustering(fi, pr, reps):
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

    mat = numpy.zeros((len(contigs),len(contigs)), dtype=numpy.float)

    for i in range(reps):
        os.system("perl makeRandom.pl {!s} {!s} ClusteringTest/dbIn.fna ClusteringTest/seqs.fna".format(pr, fi))
        os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 2 -i ClusteringTest/dbIn.fna -d ClusteringTest/db")
        os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 0 -i ClusteringTest/seqs.fna -d ClusteringTest/db -o ClusteringTest/res")
        os.system("rm ClusteringTest/dbIn.fna")
        os.system("rm ClusteringTest/seqs.fna")
        os.system("rm ClusteringTest/db")
        
        res = open("ClusteringTest/res")
        lines = res.readlines()
        for sq, db in zip(lines[1::2], lines[2::2]):
            a = sq.rstrip()[1:]
            b = db.rstrip()
            mat[contigs.index(a)][contigs.index(b)] += 1.0
        res.close()
    os.system("rm ClusteringTest/res")

    mat = mat/reps # Normalizing

    adj = open("ClusteringTest/adjacency_{!s}.csv".format(reps),'w')
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

    gra = open("ClusteringTest/graph_facts_{!s}.txt".format(reps),'w')
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
    plt.savefig("ClusteringTest/graph_{!s}.png".format(reps))
    gra.close()

def adjMatrix2Graph(fi):
    import networkx as nx
    
    G = nx.Graph()
    file = open(fi,'r')
    lines = file.readlines()
    for x in range(len(lines)):
        li = lines[x].rstrip().split(",")
        for y in range(x):
            w = float(li[y])
            if w > 0.0:
                G.add_edge(x,y,weight=float(li[y]))
    return G

def adjMatrix2Dot(num, fi, out):
    dot = open(out,'w')
    dot.write("digraph G {\n")
    file = open(fi,'r')
    lines = file.readlines()
    file.close()
    for x in range(len(lines)):
        li = lines[x].rstrip().split(",")
        for y in range(len(lines)):
            w = float(li[y])
            if w > 0.0:
                dot.write("\t{!s} -> {!s} [weight={!s}]\n".format(x,y,int(num*w)))
    dot.write("}")
    dot.close()

def Fasta2Tab(fIn,fOut):
    fileIn = open(fIn,'r')
    fileOut = open(fOut,'w')
    buf = fileIn.readline().rstrip()
    while buf:
        if buf[0] == '>':
            fileOut.write(buf)
        else:
            fileOut.write("\t{!s}\n".format(buf))
        buf = fileIn.readline().rstrip()
    fileIn.close()
    fileOut.close()

def bootStrap(fi, threshold):
    import os
    import numpy
    import scipy

    #contigs = []
    #input = open(fi,'r')
    #buf = input.readline()
    #while buf:
    #    contigs.append(buf.split('\t')[0])
    #    buf = input.readline()
    #input.close()

    #mat = numpy.zeros((len(contigs),len(contigs)), dtype=numpy.float)

    os.system("perl separateBySize.pl {!s} {!s} BootstrappingTest/dbIn.fna BootstrappingTest/seqs.fna BootstrappingTest/nums.txt".format(threshold, fi))
    nmf = open("BootstrappingTest/nums.txt",'r')
    [dbNum, seqNum] = nmf.readline().rstrip().split(",")
    os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 2 -i BootstrappingTest/dbIn.fna -d BootstrappingTest/db")
    os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 0 -i BootstrappingTest/seqs.fna -d BootstrappingTest/db -o BootstrappingTest/resSeeds")
    os.system("rm BootstrappingTest/dbIn.fna")
    os.system("rm BootstrappingTest/seqs.fna")
    os.system("rm BootstrappingTest/db")
    os.system("rm BootstrappingTest/nums.txt")

    seqs = []
    matches = []
    res = open("BootstrappingTest/resSeeds",'r')
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

    fout = open("BootstrappingTest/Clusters_{!s}.txt".format(threshold),'w')
    fout.write(str(list1))
    fout.close()


def adjMatrix2Dendrograms(fi, namefi, pathout, num):
    import numpy
    import scipy
    import scipy.cluster.hierarchy as sch
    import matplotlib.pyplot as plt
    import pickle
    import math
    import pylab
    
    if pathout[-1] != "/":
        pathout = pathout + "/"
    
    pathout = pathout + str(num)
    
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
    fig.savefig("{!s}_dendrogram_matrix.png".format(pathout))

    order = open("{!s}_dendro_order.p".format(pathout),'wb')
    pickle.dump(Z1['leaves'],order)
    order.close()

    clustersOut = open("{!s}_dendro_clusters.txt".format(pathout),'w')

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
    import math
    
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
            
        print clDict
        print clDict
                    
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


def distanceDistributionCSV(file, out):
    # http://stackoverflow.com/a/5328669
    import matplotlib.pyplot as plt
    import numpy as np
    
    f = open(file,'r')
    lines = f.readlines()
    list = []
    
    for i in range(len(lines)):
        row = lines[i].rstrip().split(',')
        for j in range(len(lines)):
            if float(row[j]) > 0:
                list.append(float(row[j]))
    f.close()

    mat = np.array(list)

    hist, bins = np.histogram(mat, bins=50)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align='center', width=width)
    plt.show()


def clusterSimilarityK(A, B):
    import numpy
    import math
    '''From "A Method for Comparing Two Hierarchical Clusterings", E.B. Fowlkes and C.L. Mallows, 1983
        Requires that A, B have the same number of clusters and the same number of total objects'''
    # Determine M where m_ij = the number of objects between the ith cluster of A and the jth cluster of B
    k = len(A)
    M = numpy.zeros((k, k), dtype=numpy.float)
    n = 0
    for i in range(k):
        cA = set(A[i])
        for j in range(k):
            cB = set(B[j])
            M[i][j] = len(cA.intersection(cB))
        n += len(A[i])

    # Determine Tk, Pk, Qk
    n = float(n)

    Mi_ = numpy.sum(M, axis=1)
    M_j = numpy.sum(M, axis=0)

    Tk = numpy.sum(M**2)-n
    Pk = numpy.sum(Mi_**2)-n
    Qk = numpy.sum(M_j**2)-n

    Pkp, Qkp = 0.0, 0.0
    for l in range(k):
        mi = Mi_[l]
        mj = M_j[l]
        Pkp += mi * (mi-1.0) * (mi-2.0)
        Qkp += mj * (mj-1.0) * (mj-2.0)

    Bk = Tk / math.sqrt(Pk*Qk)
    EBk = math.sqrt(Pk*Qk)/(n*(n-1.0))
    var = 2.0/(n*(n-1.0)) + (4.0*Pkp*Qkp)/(n*(n-1.0)*(n-2.0)*Pk*Qk) + (Pk - 2.0 - 4.0*(Pkp/Pk))*(Qk - 2.0 -4.0*(Qkp/Qk))/(n*(n-1.0)*(n-2.0)*(n-3.0)) - (Pk * Qk)/((n**2)*((n-1.0)**2))

    return Bk, EBk, math.pow(var,0.5), [EBk - 2*(math.pow(var,0.5)), EBk + 2*(math.pow(var,0.5))]

def csv2NumpyMatrix(fi):
    import numpy

    f = open(fi,'r')
    lines = f.readlines()
    mat = numpy.zeros((len(lines),len(lines)), dtype=numpy.float)
    
    for i in range(len(lines)):
        row = lines[i].rstrip().split(',')
        for j in range(len(lines)):
            mat[i][j] = float(row[j])
    f.close()

    return mat

def adjMatrix2SymMatrix(fi):
    mat = csv2NumpyMatrix(fi)
    for i in range(len(mat)):
        for j in range(i):
            m = (mat[i][j] + mat[j][i])/2
            mat[i][j] = m
            mat[j][i] = m
    return mat


def adjMatrix2SymDistanceMatrix(A):
    import numpy
    '''Takes a numpy adjacency matrix A and returns a symmetric distance matrix B'''

    M = numpy.zeros((len(A),len(A)),dtype=numpy.float)
    max = 0.0
    for i in range(len(A)):
        for j in range(i+1):
            mij = (A[i][j] + A[j][i])/2.0
            M[i][j] = mij
            M[j][i] = mij
            if mij > max:
                max = mij

    O = numpy.ones((len(A),len(A)),dtype=numpy.float)
    O = max * O

    return O - M

def getDistances(fi, pctrep):
    import numpy
    import math
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
    print numrep
            
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

    fig = plt.figure()
    n, bins, patches = plt.hist(dist, 50, facecolor='green')

    plt.savefig("{!s}Histogram.png".format(folder))

def testSelfRecruitment(pathIn, numSpecies, sizeChop, numDB, numSeq):
    import os
    import random
    import numpy as np
    import scipy.cluster.vq as scv
    import scipy.cluster.hierarchy as sch
    import networkx as nx
    import community as cm
    import matplotlib.pyplot as plt

    if pathIn[-1] != '/':
        pathIn = pathIn + '/'
    allList = os.listdir(pathIn)
    genomeList = []
    for file in allList:
        if file[-4:] == ".fna":
            genomeList.append(file)

    subset = random.sample(genomeList, numSpecies)

    for file in subset:
        # make sequences to be matched
        ensureDir(pathIn+"Sequences/")
        chopRandom(pathIn+file, pathIn+"Sequences/", sizeChop, sizeChop/10, numSeq)
        # make sequences to be matched to
        ensureDir(pathIn+"DataBase/")
        chopRandom(pathIn+file, pathIn+"DataBase/", sizeChop, sizeChop/10, numDB)

    # Make RAI databases
    os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 2 -I {!s}Sequences/ -d {!s}seqs".format(pathIn,pathIn))        
    os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 2 -I {!s}DataBase/ -d {!s}db".format(pathIn,pathIn))

    # Data sets for further evaluation
    namesSq, RaiSq = rai2Numpy(pathIn+"seqs")
    namesDb, RaiDb = rai2Numpy(pathIn+"db")
            
    namesAll = namesSq + namesDb
    RaiAll = np.concatenate([RaiSq,RaiDb])

    csvout = open("{!s}seqs.csv".format(pathIn),'w')
    for r in RaiAll:
        csvout.write(",".join(str(x) for x in r)+"\n")
            
    # Run RAIphy
    os.system("/Users/anna/Research/RAIphyCommandLine/raiphy -e .fna -m 0 -I {!s}Sequences/ -d {!s}db -o {!s}output".format(pathIn,pathIn,pathIn))
    # Evaluate RAIphy results
    raiDict = {}
    for k in namesDb:
        raiDict[k] = [k]
    raiRes = open("{!s}output".format(pathIn),'r')
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

    os.system("rm -r {!s}Sequences/".format(pathIn))
    os.system("rm -r {!s}seqs".format(pathIn))
    os.system("rm -r {!s}DataBase/".format(pathIn))
    os.system("rm -r {!s}db".format(pathIn))
    os.system("rm -r {!s}output".format(pathIn))

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
    mylen = len(namesDb)
    clustWeightedLink = sch.fcluster(WeightedLink,numSpecies,criterion='maxclust')
    hist, bins = np.histogram(clustWeightedLink, bins=numSpecies)
    print hist

# http://stackoverflow.com/a/17622583
def newWhiten(obs):
    import numpy
    std_dev = numpy.std(obs)
    return obs/std_dev

def rai2Numpy(fi):
    import numpy
    names = []
    obs = []
    
    f = open(fi,'r')
    f.readline()
    buf = f.readline().rstrip()
    while buf:
        l = buf.split(":")
        names.append(l[1][1:])
        o = [0.0]*(len(l)-2)
        for i in range(len(l)-2):
            o[i] = float(l[i+2])
        obs.append(o)
        
        buf = f.readline().rstrip()
    
    return names, numpy.array(obs)