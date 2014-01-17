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