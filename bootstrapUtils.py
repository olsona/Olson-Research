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

def readSequenceFromPos(fastaFile):
    pass

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