from decimal import *

#----For use in Bootstrap----#
def checkCorrectMatchClusterMaxSpecies(match, cluster, names):
    import string
    _, name = purityOfCluster(cluster.getAll(), names)
    if name == '-':
   	return 2
    elif string.find(match, name) != 1:
	return 1
    else:
	return 0
		
def checkCorrectMatchClusterMaxGenus(match, cluster, names):
    import string
    _, name = purityOfCluster(cluster.getAll(), names)
    if name == '-':
   	return 2
    elif string.split(match,'_')[0]==string.split(name,'_')[0]:
	return 1
    else:
	return 0
		
		
def checkCorrectMatchKnown(match, name):
    import string
    if string.find(match,name) != 1:
	return 1
    else:
	return 0
		
		
def checkCorrectSpeciesOlsonFormat(match,contig):
    import string
    if string.rsplit(match,'_',2)[0]==string.rsplit(contig,'_',2)[0]:
	return 1
    else:
	return 0
		
		
def checkCorrectGenusOlsonFormat(match,contig):
    import string
    if string.split(match,'_')[0]==string.split(contig,'_')[0]:
	return 1
    else:
	return 0
		
		
def comparisonPlot(rdata, wdata, iterString, outputFile, suffix, rlabel, wlabel):
    import matplotlib.pyplot as plt
    import numpy
	
    if rdata and wdata:
	mymax = max(max(rdata),max(wdata))
	mymin = min(min(rdata),min(wdata))
    elif rdata:
   	mymax = max(rdata)
	mymin = min(rdata)
    elif wdata:
	mymax = max(wdata)
	mymin = min(wdata)
    else:
	return
		
    bins = numpy.linspace(mymin,mymax,40)

    if rdata:
	plt.hist(rdata,bins, normed=0, facecolor='blue', alpha=0.5, label=rlabel)
    if wdata:
	plt.hist(wdata,bins, normed=0, facecolor='red', alpha=0.5, label=wlabel)
    plt.xlabel("Score")
    plt.title("Correct vs Incorrect Scores UNNORMED, {!s}".format(iterString))
    plt.legend()
    plt.savefig("{!s}_{!s}_{!s}_unnorm.pdf".format(outputFile, iterString, suffix), bbox_inches='tight')
    plt.clf()		
		
		
#----Evaluation post-computation----#
def purityOfCluster(clusterElements, nameList):
    '''Returns the ratio representation and identity of most common name out of the whole cluster, where the possible names are given in nameList.'''
    import string
    repDict = {nL: 0 for nL in nameList}
    max = 0
    maxName = ''
    for c in clusterElements:
	for nL in nameList:
	    if string.find(c,nL) != -1:
	        repDict[nL] += 1
    for nL in repDict:
	if repDict[nL] > max:
	    max = repDict[nL]
	    maxName = nL
    ratio = float(max)/float(len(clusterElements))
    return ratio, maxName
	
def purityOfClusterByLength(cluster, nameList):
    '''Returns the ratio representation (by length) and identity of most common name out of the whole cluster, where the possible names are given in nameList.'''
    if cluster is None:
	return 0.0, '-'
    else:
	import string
	repDict = {nL: 0 for nL in nameList}
	maxLen = 0
	totalLen = 0
	maxName = ''
	for c in cluster:
	    len = int(c.split("_")[-2])
	    totalLen += len
	    for nL in nameList:
	        if string.find(c,nL) != -1:
	            repDict[nL] += len
	for nL in repDict:
	    if repDict[nL] > maxLen:
	        maxLen = repDict[nL]			
	        maxName = nL
	ratio = float(maxLen)/float(totalLen)
	return ratio, maxName
	
	
def sensitivity(inFile, correctClustering, threshold, nameList):
    # tacoa_review, eqn (8)
    import string
    f = open(inFile, 'r')
    lns = f.readlines()
    ZDict = {name:0 for name in nameList}
    TPDict = {name:0 for name in nameList}
    # get number of contigs that *should* be classified for each name (Zi)
    for cl in correctClustering:
	rep = cl[0]
	corName = ''
	for name in nameList:
	    if string.find(rep,name) != -1:
	        corName = name
	        break
	Z[corName] = len(cl)
    # iterate through inFile, find clusters that have purity above a threshold, and count true positives (TPi)
    for li in lns:
        repDict = {nL: 0 for nL in nameList}
	max = 0
	maxName = ''
	fstBrk = li.rstrip().split(": ")
	ls = fstBrk[1].split(", ")
	for l in ls:
	    for nL in nameList:
		if string.find(l,nL) != -1:
			repDict[nL] += 1
	for nL in repDict:
	    if repDict[nL] > max:
		max = repDict[nL]
		maxName = nL
	ratio = float(max)/float(len(ls))
	if ratio >= threshold:
	    TPDict[maxName] += max
	# compute sensitivity
	SnDict = {name: 0.0 for name in nameList}
	for name in nameList:
	   SnDict[name] = TPDict[name]/ZDict[name]
	return SnDict
	
	
def precision(inFile, threshold, nameList):
    # tacoa_review, eqn (8)
    import string
    f = open(inFile, 'r')
    lns = f.readlines()
    TPDict = {name:0 for name in nameList}
    PDict = {name:0 for name in nameList}
    # iterate through inFile, find clusters that have purity above a threshold, and count true positives (TPi)
    for li in lns:
        repDict = {nL: 0 for nL in nameList}
        max = 0
	maxName = ''
	fstBrk = li.rstrip().split(": ")
	ls = fstBrk[1].split(", ")
	for l in ls:
	    for nL in nameList:
		if string.find(l,nL) != -1:
		    repDict[nL] += 1
	for nL in repDict:
	    if repDict[nL] > max:
	        max = repDict[nL]
	        maxName = nL
	ratio = float(max)/float(len(ls))
	if ratio >= threshold:
	    TPDict[maxName] += max
	    PDict[maxName] += len(ls)
    # compute sensitivity
    PrDict = {name: 0.0 for name in nameList}
    for name in nameList:
	PrDict[name] = TPDict[name]/PDict[name]
    return PrDict
	
	
def purityClusterWholeOutput(inFile, nameList, outFile):
    import string
    f = open(inFile, 'r')
    out = open(outFile, 'w')
    lns = f.readlines()
    for li in lns:
        repDict = {nL: 0 for nL in nameList}
	max = 0
	maxName = ''
	fstBrk = li.rstrip().split(": ")
	nm = fstBrk[0]
	ls = fstBrk[1].split(", ")
	for l in ls:
	    for nL in nameList:
		if string.find(l,nL) != -1:
		    repDict[nL] += 1
	for nL in repDict:
	    if repDict[nL] > max:
	        max = repDict[nL]
		maxName = nL
	ratio = float(max)/float(len(ls))
	out.write("{!s}:\t{:.2%}\t{!s}\n".format(nm, ratio, maxName))
    f.close()
    out.close()
	
	
def Entropy(U):
    import math
    # http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
    N = 0
    for u in U:
	N += len(u)
    sum = 0.0
    for u in U:
	p = float(len(u))/float(N)
	sum += p * math.log(p)
    return sum*(-1)
	
	
def MutualInformation(U, V):
    # http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
    import math
    K = len(U)
    J = len(V)
    N = 0
    for k in range(K):
	N += len(U[k])
    sum = 0.0
    for k in range(K):
	uk = len(U[k])
	for j in range(J):
	    vj = len(V[j])
	    s1 = set(U[k])
	    s2 = set(V[j])
	    ukvj = len(set.intersection(s1,s2))
	    A = float(ukvj)/float(N)
	    B = float(N*ukvj)/float(uk*vj)
	    if A != 0.0 and B != 0.0:
		sum += A*math.log(B)
    return sum
	
	
def NMI(U, V):
    # http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
    I = MutualInformation(U, V)
    HU = Entropy(U)
    HV = Entropy(V)
    return (2.0*I)/(HU+HV)
	
	
def ExpectedMutualInformation(U,V):
    # Vinh, Epps, Bailey, (2)
    import math
    getcontext().prec = 400
    N = 0
    for u in U:
        N += len(u)
    R = len(U)
    C = len(V)
    E = Decimal(0)
    for i in range(R):
	for j in range(C):
	    ai = len(U[i])			
	    bj = len(V[j])
	    for nij in range(max(ai+bj-N,1), min(ai,bj)+1):
	        if ai != 0 and bj != 0:
	            t1 = Decimal((float(nij)/float(N)) * (math.log(float(N*nij)/float(ai*bj))))
	            t2 = Decimal(ncr(N,nij))
	            t2 /= Decimal(ncr(N,ai))
	            t2 *= Decimal(ncr(N-nij,ai-nij))
	            t2 /= Decimal(ncr(N,bj))
	            t2 *= Decimal(ncr(N-ai,bj-nij))
	            E += t1*t2	
    return E
	
	
def ncr(n, r):
    # http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
    import operator as op
    r = min(r, n-r)
    if r == 0: return 1
    numer = reduce(op.mul, xrange(n, n-r, -1))
    denom = reduce(op.mul, xrange(1, r+1))
    return numer//denom
	
	
def AMI(U, V):
    # Information Theoretic Measures for Clusterings, Vinh, Epps, Bailey
    # AMI_sum
    I = MutualInformation(U, V)
    HU = Entropy(U)
    HV = Entropy(V)
    E = ExpectedMutualInformation(U,V)
    return (I-E)/(0.5*(HU+HV)-E)
	
	
#----Total post-processing computation----#
def testCorrectnessAll(computedClustering, correctClustering, names, outFile, representationThreshold=0.9):
    purityInfo = {c:[] for c in computedClustering}
    totalR = 0
    totalN = 0
    # check purity for each cluster
    for c in computedClustering:
	#print c
	#print computedClustering[c]
	pur, max = purityOfCluster(computedClustering[c], names)
	purityInfo[c] = [pur, max] 
	l = len(computedClustering[c])
	totalR += int(pur*l+0.5)
	totalN += l
		
    # total purity
    totalPurity = float(totalR)/float(totalN)
	
    # check representation
    repDict = {n:0 for n in names}
    for c in purityInfo:
	[pur, max] = purityInfo[c]
	if pur >= representationThreshold:
	    repDict[max] += 1
		
    # check NMI
    nmi = NMI(computedClustering.values(), correctClustering)
	
    outF = open(outFile, 'w')
    outF.write("NMI:\t{0:.6f}\n".format(nmi))
    outF.write("Total purity:\t{:03.2f}%\n".format(totalPurity*100.0))
    for c in purityInfo:
	outF.write("Purity of {!s}:\t{:03.2f}%, {!s}\n".format(c,purityInfo[c][0]*100.0,purityInfo[c][1]))
    for n in repDict:
	outF.write("Representation of {!s}:\t{!s}\n".format(n,repDict[n]))
		
    outF.close()