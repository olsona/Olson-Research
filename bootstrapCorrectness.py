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


def Entropy(clustering):
	# http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
	N = 0
	for c in clustering:
		N += len(c)
	sum = 0.0
	for c in clustering:
		p = float(len(c))/float(N)
		sum += p * log(p)
	return sum*(-1)


def MututalInformation(clustering, classes):
	# Ground truth is stored in 'classes'
	# http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
	K = len(clustering)
	J = len(classes)
	N = 0
	for k in range(K):
		N += len(clustering[k])
	sum = 0.0
	for k in range(K):
		wk = clustering[k]
		for j in range(J):
			cj = classes[j]
			wkcj = len(set.intersection(wk,cj))
			A = float(len(wkcj))/float(N)
			B = float(N*len(wkcj))/float(len(wk)*len(cj))
			sum += A*log(B)
	return sum



def NMI(clustering, classes):
	# http://nlp.stanford.edu/IR-book/html/htmledition/evaluation-of-clustering-1.html
	


def RandIndex(correct, computed):
	pass
	#return (truePos+trueNeg)/(truePos+trueNeg+falsePos+falseNeg)


def NMI(correct, computed):
	pass