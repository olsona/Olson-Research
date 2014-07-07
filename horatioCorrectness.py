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
	# aka specificity
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
	#aka specificity
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
	

def sensitivityCluster(inCluster, correctClustering, threshold, nameList):
	# tacoa_review, eqn (8)
	import string
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
	for cl in inCluster:
		repDict = {nL: 0 for nL in nameList}
		for c in cl:
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
	
	
def sensitivityPrintedFile(inFile, correctClustering, threshold, nameList):
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
	#aka specificity
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
	getcontext().prec = 20
	#getcontext().Emax = 2000
	N = 0
	for u in U:
		N += len(u)
	R = len(U)
	print R
	C = len(V)
	print C
	E = Decimal(0.0)
	start = time.time()
	for i in range(R):
		if i % 10 == 0:
			print "	 i=", i, " runtime=", (time.time() - start)
			start = time.time()
	for j in range(C):
		ai = len(U[i])
		if ai:
			#print "ai: ", N, ai;
			ncr_Nai = Decimal(ncr(N,ai))
			#print "{:2.2e}".format(ncr_Nai)			
		bj = len(V[j])
		if bj:
			#print "bj:", N, bj
			ncr_Nbj = Decimal(ncr(N,bj))
			#print "{:2.2e}".format(ncr_Nbj)
		for nij in range(max(ai+bj-N,1), min(ai,bj)+1):
			if ai != 0 and bj != 0:
				t1 = Decimal((float(nij)/float(N)) * (math.log(float(N*nij)/float(ai*bj))))
				t2 = Decimal(ncr(N,nij))
				t2 /= ncr_Nai
				t2 *= Decimal(ncr(N-nij,ai-nij))
				t2 /= ncr_Nbj
				t2 *= Decimal(ncr(N-ai,bj-nij))
				E += t1*t2	

	return E
	
ncr_memo = {}	
def ncr(n, r):
	# http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
	import operator as op
	#return 1.0
	#print n, r
	r = min(r, n-r)
	if r == 0: return 1
		nrtuple = tuple([n,r])
	if nrtuple not in ncr_memo:
		numer = reduce(op.mul, xrange(n, n-r, -1))
		denom = reduce(op.mul, xrange(1, r+1))
		ncr_memo[nrtuple] = numer//denom
	#print "{:2.2e}".format(ncr_memo[nrtuple])
	return ncr_memo[nrtuple]
	#return 1.0
	
	
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
	
	
# to process an entire folder	 
def processFolder(inFolder, nameFile, correctFilePrefix, sizeThreshold, outFile, maxsize):
	import glob, string, time, numpy
	import cPickle as pickle
	from horatioClasses import Cluster
	from scipy.stats import pearsonr
	
	cutDict = {'allClose': [2,4,6,8,10,12,14,16,18,20,22],
			'lowClose': [2,4,6,8,10,14,18,22],
			'2by4': [2,6,10,14,18,22],
			'4by4': [4,8,12,16,20],
			'4allClose': [4,6,8,10,12,14,16,18,20,22],
			'4lowClose': [4,6,8,10,14,18,22],
			'allCloseChop': [2,4,6,8,10,12,14,16,18],
			'lowCloseChop': [2,4,6,8,10,14,18],
			'2by4Chop': [2,6,10,14,18],
			'4allCloseChop': [4,6,8,10,12,14,16,18],
			'4lowCloseChop': [4,6,8,10,14,18],
			'4by4Chop': [4,8,12,16]
			}
			
	names = []
	nf = open(nameFile,'r')
	for li in nf.readlines():
		nm = li.rstrip()
		if len(nm) > 0:
			names.append(nm)
			
	outF = open(outFile,'w')
	#outF.write("Source;Abundance;Score;Cut;N;J;L;SizeThreshold;AvgClustSize;MinClustSize;MaxClustSize;NMI;SnNo;SpNo;RepNo;SnLen;SpLen;RepLen;Sn4k;Sp4k;Sn6k;Sp6k;Sn8k;Sp8k;Sn10k;Sp10k\n")
	
	# get Z values for each correct clustering
	corrList = glob.glob("{!s}*".format(correctFilePrefix))
	corrDict = {int(cf.split("_gt")[1][0:-1]):cf for cf in corrList}
	corrZNo = {no:{na: 0 for na in names} for no in corrDict}
	corrZLen = {no:{na: 0 for na in names} for no in corrDict}
	Zk = {i:0 for i in [4,6,8,10]}
	for no in corrDict:
		clustering = pickle.load(open(corrDict[no],'rb'))
		for cl in clustering:
			rep = cl[0]
			corName = ''
			for nm in names:
				if string.find(rep,nm) != -1:
					corName = nm
					break
			corrZNo[no][corName] = len(cl)
			for c in cl:
				clen = int(c.rsplit("_",2)[1])
				#if clen >= 4000:
				#	 Zk[4] += 1
				#if clen >= 6000:
				#	 Zk[6] += 1
				#if clen >= 8000:
				#	 Zk[8] += 1
				#if clen >= 10000:
				#	 Zk[10] += 1
				corrZLen[no][corName] += clen
	
	ctr = 0
	threshold = 0.0
	
	fileList = glob.glob("{!s}/*_pickle".format(inFolder))
	start = time.time()
	for fi in fileList:
		clustMemList = []
		#avgScore = []
		#lowScore = []
		#sdScore = []
		purityNo = []
		purityLen = []
		#start = time.time()
		if ctr % 10 == 0:
			ntime = time.time()
			print "	  i = {:d}, {:.2f}s".format(ctr,ntime-start)
			start = time.time()
		# get information from filename
		fileName = fi.split("/")[-1]
		fileSplit = fileName.split("_")
		mText = fileSplit[0]
		mAbund = fileSplit[1]
		score = fileSplit[2]
		nInd = fileSplit.index("N")
		n = float(fileSplit[nInd+1])
		jInd = fileSplit.index("J")
		j = float(fileSplit[jInd+1])
		cInd = fileSplit.index("C")
		cText = fileSplit[cInd+1]
		cut = cutDict[cText]
		lInd = fileSplit.index("L")
		l = float(fileSplit[lInd+1])
		
		corrClustName = corrDict[cut[0]]
		corrClust = pickle.load(open(corrClustName,'rb'))
		ZDictNo = corrZNo[cut[0]]
		ZDictLen = corrZLen[cut[0]]
		
		repDictNo = {na:0 for na in names}
		repDictLen = {na:0 for na in names}
		
		inClustPre = pickle.load(open(fi,"rb"))
		inClust = []
		list2Clusters = {}
		for i in inClustPre:
			ncl = inClustPre[i].getAllLeaves()
			list2Clusters[ncl[0]] = inClustPre[i]
			inClust.append(ncl)
		#print inClust
		TPDictNo = {na:0 for na in names}
		FPDictNo = {na:0 for na in names}
		TPDictLen = {na:0 for na in names}
		FPDictLen = {na:0 for na in names}
		TPk = {i:0 for i in [4,6,8,10]}
		FPk = {i:0 for i in [4,6,8,10]}
		for cl in inClust:
			clustMemList.append(len(cl))
			myRepDictNo = {nL: 0 for nL in names}
			myRepDictLen = {nL: 0 for nL in names}
			maxNo = 0
			maxNoName = ''
			maxLen = 0
			maxLenName = ''
			totalLen = 0
			# get information about scores
			#scores = list2Clusters[cl[0]].getLeafScores()
			#if scores and len(cl) >= sizeThreshold:
			#	 avgScore.append(float(numpy.mean(scores)))
			#	 lowScore.append(float(numpy.min(scores)))
			#	 sdScore.append(float(numpy.std(scores)))
			#elif len(cl) >= sizeThreshold:
			#	 avgScore.append(0.0)
			#	 lowScore.append(0.0)
			#	 sdScore.append(0.0)
			totalK = {i:0 for i in [4,6,8,10]}
			repK = {na:{i:0 for i in [4,6,8,10]} for na in names}
			
			#get representations for each name
			for c in cl:
				#print c
				mylen = int(c.rsplit('_',2)[1])
				if mylen < maxsize:
					for nL in names:
						if string.find(c,nL) != -1:
							myRepDictNo[nL] += 1
							myRepDictLen[nL] += mylen
							totalLen += mylen
							for i in [4,6,8,10]:
								if mylen > i*1000:
									totalK[i] += 1
									repK[nL][i] += 1
							break
						
			#pprint.pprint(repDict)
			for nL in myRepDictNo:
				if myRepDictNo[nL] > maxNo:
					maxNo = myRepDictNo[nL]
					maxNoName = nL
			if len(cl) >= sizeThreshold:
				for nL in myRepDictLen:
					if myRepDictLen[nL] > maxLen:
						maxLen = myRepDictLen[nL]
						maxLenName = nL
				tpn = maxNo
				fpn = len(cl) - tpn
				tpl = maxLen
				fpl = totalLen - tpl
				pN = float(tpn)/float(len(cl))
				pL = float(tpl)/float(totalLen)
				purityNo.append(pN)
				purityLen.append(pL)
				if pN > threshold:
					repDictNo[maxNoName] += 1
					TPDictNo[maxNoName] += tpn
					FPDictNo[maxNoName] += fpn
				if pL > threshold:
					repDictLen[maxLenName] += 1
					TPDictLen[maxLenName] += tpl
					FPDictLen[maxLenName] += fpl
				for i in [4,6,8,10]:
					TPk[i] += repK[maxNoName][i]
					FPk[i] += totalK[i] - repK[maxNoName][i]
		#print
		#print
		
		#pprint.pprint(TPDict)
		
		#SnDict = {na:0 for na in names}
		#SpDict = {na:0 for na in names}
		ZAllNo = 0
		TPAllNo = 0
		FPAllNo = 0
		ZAllLen = 0
		TPAllLen = 0
		FPAllLen = 0
		for na in names:
			#if ZDict[na] and TPDict[na]:	  # tacoa p 13 "The overall specificity is computed over those classes that have a defined specificity value"
			ZAllNo += ZDictNo[na]
			TPAllNo += TPDictNo[na]
			FPAllNo += FPDictNo[na]
			ZAllLen += ZDictLen[na]
			TPAllLen += TPDictLen[na]
			FPAllLen += FPDictLen[na]
			#SnDict[na] = float(TPDict[na])/float(ZDict[na])			   # tacoa (8)
			#SpDict[na] = float(TPDict[na])/float(TPDict[na]+FPDict[na])   # tacoa (9)
		
		SnAllNo = 0.0
		SpAllNo = 0.0
		SnAllLen = 0.0
		SpAllLen = 0.0
		SnK = {i:0.0 for i in [4,6,8,10]}
		SpK = {i:0.0 for i in [4,6,8,10]}
		if TPAllNo:
			SnAllNo = float(TPAllNo)/float(ZAllNo)
			SpAllNo = float(TPAllNo)/float(TPAllNo+FPAllNo)
		if TPAllLen:
			SnAllLen = float(TPAllLen)/float(ZAllLen)
			SpAllLen = float(TPAllLen)/float(TPAllLen+FPAllLen)
		for i in [4,6,8,10]:
			if TPk[i]:
				SnK[i] = float(TPk[i])/float(Zk[i])
				SpK[i] = float(TPk[i])/float(TPk[i]+FPk[i])
		
		repNumNo = 0
		for r in repDictNo:
			if repDictNo[r]:
				repNumNo += 1
		repFracNo = float(repNumNo)/float(len(names))
		repNumLen = 0
		for r in repDictLen:
			if repDictLen[r]:
				repNumLen += 1
		repFracLen = float(repNumLen)/float(len(names))
		
		nmi = NMI(inClust, corrClust)
		#ami = AMI(inClust, corrClust)
		
		avgClustSize = float(sum(clustMemList))/float(len(clustMemList))
		minClustSize = min(clustMemList)
		maxClustSize = max(clustMemList)
		#try:
		#	 avg2pN, _ = pearsonr(avgScore,purityNo)
		#	 avg2pL, _ = pearsonr(avgScore,purityLen)
		#	 low2pN, _ = pearsonr(lowScore,purityNo)
		#	 low2pL, _ = pearsonr(lowScore,purityLen)
		#	 sd2pN, _ = pearsonr(sdScore,purityNo)
		#	 sd2pL, _ = pearsonr(sdScore,purityLen)
		#	 pN2pL, _ = pearsonr(purityNo,purityLen)
		#	 mpN = numpy.mean(purityNo)
		#	 mpL = numpy.mean(purityLen)
		#except:
		#	 print len(avgScore), len(lowScore), len(sdScore), len(purityScore), len(purityLen)
		#	 print("{!s};{!s};{!s};{!s};{:01.2f};{:01.2f};{:01.2f};".format(mText,mAbund,score,cText,n,j,l))
		
		outF.write("{!s};{!s};{!s};{!s};{:01.2f};{:01.2f};{:01.2f};".format(mText,mAbund,score,cText,n,j,l))
		outF.write("{!s};{:03.2f};{!s};{!s};".format(sizeThreshold,avgClustSize,minClustSize,maxClustSize))
		outF.write("{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}".format(nmi,SnAllNo,SpAllNo,repFracNo,SnAllLen,SpAllLen,repFracLen))
		#outF.write(";{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}\n".format(avg2pN,avg2pL,low2pN,low2pL,sd2pN,sd2pL,pN2pL,mpN,mpL))
		for i in [4,6,8,10]:
			#outF.write(";{:01.4f};{:01.4f}".format(SnK[i],SpK[i]))
			outF.write(";{:01.4f}".format(SpK[i]))
		outF.write("\n")
		ctr += 1
		
		
	outF.close()
	print "done"
	

def processDistLog(inFolder, out):
	import glob
	outF = open(out,'w')
	fileList = glob.glob("{!s}/*_distLog".format(inFolder))
	for fi in fileList:
		fileName = fi.split("/")[-1]
		fileSplit = fileName.split("_")
		mText = fileSplit[0]
		mAbund = fileSplit[1]
		score = fileSplit[2]
		nInd = fileSplit.index("N")
		n = float(fileSplit[nInd+1])
		jInd = fileSplit.index("J")
		j = float(fileSplit[jInd+1])
		cInd = fileSplit.index("C")
		cText = fileSplit[cInd+1]
		lInd = fileSplit.index("L")
		l = float(fileSplit[lInd+1])
		
		heading = "Source;Abundance;Score;Cut;N;J;L;Range;RightNum;RightMin;RightMax;WrongNum;WrongMin;WrongMax"
		
		f = open(fi,'r')
		lines = f.readlines()
		for li in lines:
			li = li.rstrip()
			if li:
				if li[0].isdigit():
					outF.write("{!s};{!s};{!s};{!s};{!s};{!s};{!s};{!s};".\
						format(mText, mAbund, score, n, j, cText, l, li.rstrip()))
				elif li[0] == "R":
					num = li.split("(")[1].split(")")[0]
					rs = li.split(" ")
					outF.write("{!s};{!s};{!s};".format(num,rs[-1],rs[-3]))
				elif li[0] == "W":
					num = li.split("(")[1].split(")")[0]
					rs = li.split(" ")
					outF.write("{!s};{!s};{!s}\n".format(num,rs[-1],rs[-3]))

	outF.close()
	f.close()


def processOne(inFile, nameFile, correctFile, sizeThreshold, outFile, maxsize):
	import string
	import cPickle as pickle
	from horatioClasses import Cluster
			
	names = []
	nf = open(nameFile,'r')
	for li in nf.readlines():
		nm = li.rstrip()
		if len(nm) > 0:
			names.append(nm)
			
	outF = open(outFile,'w')
	#outF.write("Source;Abundance;Score;Cut;N;J;L;SizeThreshold;AvgClustSize;MinClustSize;MaxClustSize;NMI;SnNo;SpNo;RepNo;SnLen;SpLen;RepLen;Sn4k;Sp4k;Sn6k;Sp6k;Sn8k;Sp8k;Sn10k;Sp10k\n")
	
	# get Z values for each correct clustering
	corrZNo = {na:0 for na in names}
	corrZLen = {na:0 for na in names}
	corrClust = pickle.load(open(correctFile,'rb'))
	for cl in corrClust:
		rep = cl[0]
		corName = ''
		for nm in names:
			if string.find(rep,nm) != -1:
				corName = nm
				break
		corrZNo[corName] = len(cl)
		for c in cl:
			clen = int(c.rsplit("_",2)[1])
			corrZLen[corName] += clen
	
	threshold = 0.0
	
	clustMemList = []
	purityNo = []
	purityLen = []
		
	repDictNo = {na:0 for na in names}
	repDictLen = {na:0 for na in names}
		
	inClustPre = pickle.load(open(inFile,"rb"))
	inClust = []
	list2Clusters = {}
	for i in inClustPre:
		ncl = inClustPre[i].getAllLeaves()
		list2Clusters[ncl[0]] = inClustPre[i]
		inClust.append(ncl)
	#print inClust
	TPDictNo = {na:0 for na in names}
	FPDictNo = {na:0 for na in names}
	TPDictLen = {na:0 for na in names}
	FPDictLen = {na:0 for na in names}
	TPk = {i:0 for i in [4,6,8,10]}
	FPk = {i:0 for i in [4,6,8,10]}
	for cl in inClust:
		clustMemList.append(len(cl))
		myRepDictNo = {nL: 0 for nL in names}
		myRepDictLen = {nL: 0 for nL in names}
		maxNo = 0
		maxNoName = ''
		maxLen = 0
		maxLenName = ''
		totalLen = 0
		totalK = {i:0 for i in [4,6,8,10]}
		repK = {na:{i:0 for i in [4,6,8,10]} for na in names}
		
		print cl
			
		#get representations for each name
		for c in cl:
			try:
				mylen = int(c.rsplit('_',2)[1])
				#print mylen
				#if mylen <= maxsize:
				for nL in names:
					if string.find(c,nL) != -1:
						myRepDictNo[nL] += 1
						myRepDictLen[nL] += mylen
						totalLen += mylen
						for i in [4,6,8,10]:
							if mylen > i*1000:
								totalK[i] += 1
								repK[nL][i] += 1
						break
			except:
				pass
		   
		if totalLen:			 
			#pprint.pprint(repDict)
			for nL in myRepDictNo:
				if myRepDictNo[nL] > maxNo:
					maxNo = myRepDictNo[nL]
					maxNoName = nL
			if len(cl) >= sizeThreshold:
				for nL in myRepDictLen:
					if myRepDictLen[nL] > maxLen:
						maxLen = myRepDictLen[nL]
						maxLenName = nL
				tpn = maxNo
				fpn = len(cl) - tpn
				tpl = maxLen
				fpl = totalLen - tpl
				pN = float(tpn)/float(len(cl))
				pL = float(tpl)/float(totalLen)
				purityNo.append(pN)
				purityLen.append(pL)
				if pN > threshold:
					repDictNo[maxNoName] += 1
					TPDictNo[maxNoName] += tpn
					FPDictNo[maxNoName] += fpn
				if pL > threshold:
					repDictLen[maxLenName] += 1
					TPDictLen[maxLenName] += tpl
					FPDictLen[maxLenName] += fpl
				for i in [4,6,8,10]:
					TPk[i] += repK[maxNoName][i]
					FPk[i] += totalK[i] - repK[maxNoName][i]
		#print
		#print
		
		#pprint.pprint(TPDict)
		
		#SnDict = {na:0 for na in names}
		#SpDict = {na:0 for na in names}
	ZAllNo = 0
	TPAllNo = 0
	FPAllNo = 0
	ZAllLen = 0
	TPAllLen = 0
	FPAllLen = 0
	ZDictNo = corrZNo
	ZDictLen = corrZLen
	for na in names:
		#if ZDict[na] and TPDict[na]:	  # tacoa p 13 "The overall specificity is computed over those classes that have a defined specificity value"
		ZAllNo += ZDictNo[na]
		TPAllNo += TPDictNo[na]
		FPAllNo += FPDictNo[na]
		ZAllLen += ZDictLen[na]
		TPAllLen += TPDictLen[na]
		FPAllLen += FPDictLen[na]
		#SnDict[na] = float(TPDict[na])/float(ZDict[na])			   # tacoa (8)
		#SpDict[na] = float(TPDict[na])/float(TPDict[na]+FPDict[na])   # tacoa (9)
		
	SnAllNo = 0.0
	SpAllNo = 0.0
	SnAllLen = 0.0
	SpAllLen = 0.0
	SnK = {i:0.0 for i in [4,6,8,10]}
	SpK = {i:0.0 for i in [4,6,8,10]}
	if TPAllNo:
		SnAllNo = float(TPAllNo)/float(ZAllNo)
		SpAllNo = float(TPAllNo)/float(TPAllNo+FPAllNo)
	if TPAllLen:
		SnAllLen = float(TPAllLen)/float(ZAllLen)
		SpAllLen = float(TPAllLen)/float(TPAllLen+FPAllLen)
	for i in [4,6,8,10]:
		if TPk[i]:
			SnK[i] = float(TPk[i])/float(Zk[i])
			SpK[i] = float(TPk[i])/float(TPk[i]+FPk[i])
		
	repNumNo = 0
	for r in repDictNo:
		if repDictNo[r]:
			repNumNo += 1
	repFracNo = float(repNumNo)/float(len(names))
	repNumLen = 0
	for r in repDictLen:
		if repDictLen[r]:
			repNumLen += 1
	repFracLen = float(repNumLen)/float(len(names))
		
	nmi = NMI(inClust, corrClust)
	#ami = AMI(inClust, corrClust)
		
	avgClustSize = float(sum(clustMemList))/float(len(clustMemList))
	minClustSize = min(clustMemList)
	maxClustSize = max(clustMemList)
		
	#outF.write("{!s};{!s};{!s};{!s};{:01.2f};{:01.2f};{:01.2f};".format(mText,mAbund,score,cText,n,j,l))
	outF.write("{!s};{:03.2f};{!s};{!s};".format(sizeThreshold,avgClustSize,minClustSize,maxClustSize))
	outF.write("{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}".format(nmi,SnAllNo,SpAllNo,repFracNo,SnAllLen,SpAllLen,repFracLen))
	#outF.write(";{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}\n".format(avg2pN,avg2pL,low2pN,low2pL,sd2pN,sd2pL,pN2pL,mpN,mpL))
	for i in [4,6,8,10]:
		#outF.write(";{:01.4f};{:01.4f}".format(SnK[i],SpK[i]))
		outF.write(";{:01.4f}".format(SpK[i]))
	outF.write("\n")
		
	outF.close()
	print "done"