#----For use in Bootstrap----#
def checkCorrectMatchClusterMax(match, cluster):
	import string
	_, name = purityOfCluster(cluster)
	if string.find(match, name) != 1:
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


#----Evaluation post-computation----#
def purityOfCluster(cluster, nameList):
	'''Returns the percentage representation and identify of most common name out of the whole cluster, where the possible names are given in nameList.'''
	import string
	repDict = {nL: 0 for nL in nameList}
	max = 0
	maxName = ''
	for c in cluster:
		for nL in nameList:
			if string.find(c,nL) != -1:
				repDict[nL] += 1
				if repDict[nL] > max:
					max =repDict[nL]
					maxName = nL
				break
	
	ratio = float(max)/float(len(cluster))
	return ratio, maxName


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
					if repDict[nL] > max:
						max = repDict[nL]
						maxName = nL
					break
		ratio = float(max)/float(len(ls))
		out.write("{!s}:\t{:.2%}\t{!s}\n".format(nm, ratio, maxName))
	f.close()
	out.close()


def purityAll(inFile, nameList):
	# http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
	import string
	f = open(inFile, 'r')
	lns = f.readlines()
	N = 0
	R = 0
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
					if repDict[nL] > max:
						max = repDict[nL]
						maxName = nL
					break
		R += max
		N += len(fstBrk)
	return float(R)/float(N)
	f.close()


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
	N = 0
	for u in U:
		N += len(u)
	R = len(U)
	C = len(V)
	E = 0.0
	for i in range(R):
		for j in range(C):
			ai = len(U[i])			
			bj = len(V[j])
			for nij in range(max(ai+bj-N,1), min(ai,bj)+1):
				if ai != 0 and bj != 0:
					t1 = (float(nij)/float(N)) * (math.log(float(N*nij)/float(ai*bj)))
					print ai, bj, N-ai, N-bj
					print N, nij, ai-nij, bj-nij, N-ai-bj+nij
					#upper = math.factorial(ai)*math.factorial(bj)*math.factorial(N-ai)*math.factorial(N-bj)
					#lower = math.factorial(N)*math.factorial(nij)*math.factorial(ai-nij)*math.factorial(bj-nij)*math.factorial(N-ai-bj+nij)
					#E += t1 * float(upper)/float(lower)	
					E += t1						
	return E


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
		pur, max = purityOfCluster(c, names)
		purityInfo[c] = [pur, max] 
		l = len(computedClustering[c])
		totalR += int(pur*l+0.5)
		totalN += l
		
	# total purity
	totalPurity = totalR/totalN
	
	# check representation
	repDict = {n:0 for n in names}
	for c in purityInfo:
		[pur, max] = purityInfo[c]
		if pur >= representationThreshold:
			repDict[max] += 1
		
	# check NMI
	nmi = NMI(list(computedClustering.values()), list(correctClustering.values()))
	
	outF = open(outFile, 'w')
	outF.write("NMI:\t%0.8f\n".format(nmi))
	outF.write("Total purity:\t%0.8f\n".format(totalPurity))
	for c in purityInfo:
		outF.write("Purity of %s:\t%0.8f,%s\n".format(c,purityInfo[c][0],purityInfo[c][1]))
	for n in repDict:
		outF.write("Representation of %s:\t%d\n".format(n,repDict[n]))
		
	outF.close()