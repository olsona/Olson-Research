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


def spectralContrastAngle(v):
    # From J Am Soc Mass Spectrom 2002, 13, 85-88, K. Wan, I. Vidavsky, M. Gross; eqn 5
    from numpy.linalg import norm
    from numpy import dot
    import math
    [vec1,vec2] = v
    v1n = norm(vec1)
    v2n = norm(vec2)
    c = dot(vec1,vec2)
    return math.acos(c/(v1n*v2n))


def euclidean(v):
    from numpy.linalg import norm
    return LA.norm(v[0] - v[1])


def negDist2(v):
    from numpy.linalg import norm
    from math import pow
    return -pow(norm(v[0]-v[1]),2)


def similarityIndex(v):
    # From J Am Soc Mass Spectrom 2002, 13, 85-88, K. Wan, I. Vidavsky, M. Gross; eqn 2
    import math
    [vec1,vec2] = v
    a = 0.0
    a0 = 0.0
    s = 0.0
    for j in range(len(vec1)):
        a = max(math.fabs(vec1[j]),math.fabs(vec2[j]))
        a0 = min(math.fabs(vec1[j]),math.fabs(vec2[j]))
        if a != a0:
            s += math.pow(100.0*(a-a0)/(a+a0), 2.0)
    return math.sqrt(s/len(vec1))

def raiScore(v):
    import math
    [vec1,vec2] = v
    s = 0.0
    for i in range(len(vec1)):
        s += vec1[i]*vec2[i]
    return s