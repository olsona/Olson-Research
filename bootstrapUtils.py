def csv2MediansTop(inFile, outFile, numTop):
    import numpy
    mat = numpy.genfromtxt(inFile, dtype=numpy.float32, delimiter=" ")
    fout = open(outFile, 'w')
    for r in range(len(mat)/numTop):
        segment = mat[r*numTop:(r+1)*numTop]
        fout.write("{!s}\n".format(numpy.median(segment,axis=0)[2]))
    fout.close()


def ensureDir(pth):
    # source: http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    import os
    d = os.path.dirname(pth)
    if not os.path.exists(d):
        os.makedirs(d)


def namesPosTable(fastaFile):
    import re, string
    f = open(fastaFile)
    table = {}
    ln = f.readline()
    while ln:
        if ln[0] == '>': # found a contig name
            nm = string.replace(ln.rstrip(),' ','_')
            m = re.search('[A-Za-z]',ln).start()
            table[nm[m:]] = f.tell()-len(ln)
        ln = f.readline()
    return table


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


def getLeaves(inDict, root):
    leaves = []
    for l in inDict[root]:
        if l in inDict.keys():
            res = getLeaves(inDict,l)
            for r in res:
                leaves.append(r)
        else:
            leaves.append(l)
    return leaves


def purityOfCluster(cluster, nameList):
    '''Returns the percentage representation and identify of most common name out of the whole cluster, where the possible names are given in nameList.'''
    import string
    repDict = {nL: 0.0 for nL in nameList}
    for c in cluster:
        for nL in nameList:
            if string.find(c,nL) != -1:
                repDict[nL] += 1.0
                break
    max = 0.0
    maxName = ''
    for r in repDict:
        if repDict[r] > max:
            max = repDict[r]
            maxName = r

    ratio = max/float(len(cluster))
    return ratio, maxName


def purityWholeOutput(inFile, nameList, outFile):
    import string
    f = open(inFile, 'r')
    out = open(outFile, 'w')
    lns = f.readlines()
    for li in lns:
        repDict = {nL: 0.0 for nL in nameList}
        max = 0.0
        maxName = ''
        fstBrk = li.rstrip().split(": ")
        nm = fstBrk[0]
        ls = fstBrk[1].split(", ")
        for l in ls:
            for nL in nameList:
                if string.find(l,nL) != -1:
                    repDict[nL] += 1.0
                    if repDict[nL] > max:
                        max = repDict[nL]
                        maxName = nL
                    break
        ratio = max/float(len(cluster))
        out.write("{!s}: {!f}% {!s}\n".format(nm, ratio, maxName))
    f.close()
    out.close()


def makeDistanceMatrix(scoreFile):
    import numpy
    f = open(scoreFile,'r')
    lns = f.readlines()
    rows = len(lns)-2
    cols = len(lns[2].split(", "))
    dists = numpy.zeros((rows,cols))
    for r in range(rows):
        li = lns[r+2].split(", ")
        for l in li:
            [sStr, cStr] = l.split(":")
            s = double(sStr)
            c = int(cStr)
            dists[r,c] = s
    return dists