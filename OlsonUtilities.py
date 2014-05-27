def readSequence(fi):
    '''This assumes that only one organism is in a file'''
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


def writeProperSequences(path1, path2):
    import os
    fls = os.listdir(path1)
    for f in fls:
        if f[-5:] == 'fasta':
            f2 = open(path2 + "/" + f, 'w')
            nm, se = readSequence(path1 + "/" + f)
            f2.write(">: " + nm + "\n")
            ct = 0
            while 1:
                f2.write(se[ct:ct+60] + "\n")
                ct += 60
                if (ct > len(se)):
                    break
            f2.close()


def ensureDir(pth):
    # source: http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    import os
    d = os.path.dirname(pth)
    if not os.path.exists(d):
        os.makedirs(d)


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


def csv2NumpyMatrix(fi):
    import numpy
    
    f = open(fi,'r')
    lines = f.readlines()
    mat = numpy.zeros((len(lines),len(lines)), dtype=numpy.float)
    
    for i in range(len(lines)):
        row = lines[i].rstrip().split(',')
        for j in range(len(row)):
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


def spectralContrastAngle(vec1,vec2):
    # From J Am Soc Mass Spectrom 2002, 13, 85-88, K. Wan, I. Vidavsky, M. Gross; eqn 5
    from numpy.linalg import norm
    from numpy import dot
    import math
    v1n = norm(vec1)
    v2n = norm(vec2)
    c = dot(vec1,vec2)
    return math.acos(c/(v1n*v2n))


def euclidean(vec1, vec2):
    from numpy.linalg import norm
    return norm(vec1 - vec2)


def negDist2(vec1, vec2):
    from numpy.linalg import norm
    from math import pow
    return -pow(norm(vec1-vec2),2)


def similarityIndex(vec1, vec2):
    # From J Am Soc Mass Spectrom 2002, 13, 85-88, K. Wan, I. Vidavsky, M. Gross; eqn 2
    import math
    a = 0.0
    a0 = 0.0
    s = 0.0
    for j in range(len(vec1)):
        a = max(math.fabs(vec1[j]),math.fabs(vec2[j]))
        a0 = min(math.fabs(vec1[j]),math.fabs(vec2[j]))
        if a != a0:
            s += math.pow(100.0*(a-a0)/(a+a0), 2.0)
    return math.sqrt(s/len(vec1))


def grayCode(base, digits, value):
    baseN = [0]*digits
    gray = [0]*digits
    for i in range(digits):
        baseN[i] = value % base
        value = value / base
    shift = 0
    while (i >= 0):
        gray[i] = (baseN[i] + shift) % base
        shift = shift + base - gray[i]
        i -= 1
    return gray


def csv2CoordinateInput(inFile, outFile):
    fin = open(inFile,'r')
    fout = open(outFile, 'w')
    count = 0
    buf = fin.readline().rstrip()
    while buf:
        l = buf.split(",")
        for j in range(len(l)):
            if j != count:
                fout.write("{!s} {!s} {!s}\n".format(count,j,l[j]))
        count += 1
        buf = fin.readline().rstrip()
    fin.close()
    fout.close()


def csv2MediansAll(inFile, outFile):
    import numpy
    mat = csv2NumpyMatrix(inFile)
    fout = open(outFile, 'w')
    for r in mat:
        fout.write("{!s}\n".format(numpy.median(r)))
    fout.close()


def csv2MediansTop(inFile, outFile, numTop):
    import numpy
    mat = numpy.genfromtxt(inFile, dtype=numpy.float32, delimiter=" ")
    fout = open(outFile, 'w')
    for r in range(len(mat)/numTop):
        segment = mat[r*numTop:(r+1)*numTop]
        fout.write("{!s}\n".format(numpy.median(segment,axis=0)[2]))
    fout.close()


def reverseComplement(seq):
    rules = {'a':'t', 'c':'g', 'g':'c', 't':'a'}
    op = ''
    for c in seq:
        op = rules[c]+op
    return op


def kmerCountSequence(seq, k):
    myDict = {}
    count = 0
    gc = 0
    seqL = seq.lower()
    for c in range(len(seqL)-k+1):
        kseq = seqL[c:c+k]
        if kseq in myDict:
            myDict[kseq] += 1
            myDict[reverseComplement(kseq)] += 1
        else:
            myDict[kseq] = 1
            myDict[reverseComplement(kseq)] = 1
        count += 2
        ch = seqL[c]
        if ch == 'g' or ch == 'c':
            gc += 2
    for c in range(len(seq)-k+1, len(seq)):
        ch = seqL[c]
        if ch == 'g' or ch == 'c':
            gc += 2
    return myDict, count, gc


def sizeString(n):
    if n >= 1000:
        return str(n/1000) + "k"
    else:
        return str(n)


def setBoxColors(bp):
    # http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots
    setp(bp['boxes'][0], color='blue')
    setp(bp['caps'][0], color='blue')
    setp(bp['caps'][1], color='blue')
    setp(bp['whiskers'][0], color='blue')
    setp(bp['whiskers'][1], color='blue')
    setp(bp['fliers'][0], color='blue')
    setp(bp['fliers'][1], color='blue')
    setp(bp['medians'][0], color='blue')
    
    setp(bp['boxes'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
    setp(bp['fliers'][2], color='red')
    setp(bp['fliers'][3], color='red')
    setp(bp['medians'][1], color='red')