from decimal import *
import pprint

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
		if p > 0:
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
		if uk > 0:
			for j in range(J):
				vj = len(V[j])
				if vj > 0:
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
	for i in range(R):
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
	
	
def generateLabeling(clustering):
	allItems = []
	for cl in clustering:
		for c in cl:
			allItems.append(c)
	sortItems = sorted(allItems)
	L = [0]*len(sortItems)
	ind = 0
	for cl in clustering:
		for c in cl:
			id = sortItems.index(c)
			L[id] = ind
		ind += 1
	return L	
	
	
def generateLabelingSubset(clustering, subset):
	sortItems = sorted(subset)
	L = [0]*len(sortItems)
	ind = 0
	for cl in clustering:
		for c in cl:
			if c in subset:
				id = sortItems.index(c)
				L[id] = ind
		ind += 1
	return L
	
	
# to process an entire folder	 
def processFolder(inFolder, nameFile, correctFilePrefix, sizeThreshold, outFile):
	import glob, string, time, numpy
	import cPickle as pickle
	from horatioClasses import Cluster
	from scipy.stats import pearsonr
	import sklearn.metrics
	
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
		#print li
		nm = li.rstrip()
		if len(nm) > 0:
			names.append(nm)
	#print names
	outF = open(outFile,'w')
	outF.write("Source;Score;Cut;N;J;L;Pref;SizeThreshold;NumberClusters;AvgClustSize;MinClustSize;MaxClustSize;ClustersPerGenus;NMI;AMI;V-score;SnAllNo;SpAllNo;F1No;RepFracNo;SnAllLen;SpAllLen;F1Len;repFracLen")
	#for na in sorted(names):
	#	[ge,sp,_] = na.split("_",2)
	#	name = "{0}.{1}.".format(ge[:5],sp[:5])
	#	outF.write(";{0}TPNo;{0}FPNo;{0}GrNo;{0}ZNo;{0}TPLen;{0}FPLen;{0}GrNo;{0}ZLen".format(name))
	outF.write("\n")

	# get Z values for each correct clustering
	corrList = glob.glob("{!s}*".format(correctFilePrefix))
	corrDict = {int(cf.split("gt")[1][0:-1]):cf for cf in corrList}
	corrZNo = {no:{na: 0 for na in names} for no in corrDict}
	corrZLen = {no:{na: 0 for na in names} for no in corrDict}
	for no in corrDict:
		clustering = pickle.load(open(corrDict[no],'rb'))
		#true_labels[no] = generateLabeling(clustering)
		for cl in clustering:
			if len(cl) > 0:
				rep = cl[0]
				corName = ''
				for nm in names:
					#print string.find(rep,nm), nm
					if string.find(rep,nm) != -1:
						corName = nm
						break
				#print rep, corName
				corrZNo[no][corName] = len(cl)
				#print rep, corName
				#print
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
	
	# iterate through files
	fileList = glob.glob("{!s}/*_pickle".format(inFolder))
	start = time.time()
	for fi in fileList:
		clustMemList = []
		#avgScore = []
		#lowScore = []
		#sdScore = []
		purityNo = []
		purityLen = []
		if ctr % 10 == 0:
			ntime = time.time()
			print "	  i = {:d}, {:.2f}s".format(ctr,ntime-start)
			start = time.time()
		# get information from filename
		fileName = fi.split("/")[-1]
		fileSplit = fileName.split("_")
		#print fileSplit
		mText = fileSplit[0].split(".")[0]
		#mAbund = fileSplit[1]
		score = fileSplit[1]
		nInd = fileSplit.index("N")
		n = float(fileSplit[nInd+1])
		jInd = fileSplit.index("J")
		j = float(fileSplit[jInd+1])
		cInd = fileSplit.index("C")
		cText = fileSplit[cInd+1]
		cut = cutDict[cText]
		lInd = fileSplit.index("L")
		l = float(fileSplit[lInd+1])
		pF = string.find(fi,"_A_")
		#print fi[pF:], fileSplit
		if pF >= 0:
			pInd = fileSplit.index("A")
			pref = fileSplit[pInd+1]
		else:
			pref = "N"
		if pref == "max":
			pref = "none"
		
		repDictNo = {na:0 for na in names}
		repDictLen = {na:0 for na in names}

		inClustPre = pickle.load(open(fi,"rb"))
		inClust = []
		allMembers = set()
		for i in inClustPre:
			if isinstance(i,list):
				ncl = i
			else:
				ncl = i.getLeaves()
			inClust.append(ncl)
			allMembers = allMembers | set(ncl)
			
		# because not all contigs may appear in the final clusters	
		corrClustName = corrDict[cut[0]]
		corrClustIdeal = pickle.load(open(corrClustName,'rb'))
		corrClustReal = []
		for clust in corrClustIdeal:
			clSet = set(clust)
			#clSet = set(clust) & allMembers
			corrClustReal.append(list(clSet))
		ZDictNo = corrZNo[cut[0]]
		ZDictLen = corrZLen[cut[0]]	
			
		TPDictNo = {na:0 for na in names}
		FPDictNo = {na:0 for na in names}
		TPDictLen = {na:0 for na in names}
		FPDictLen = {na:0 for na in names}
		GrDictNo = {na:0 for na in names}
		GrDictLen = {na:0 for na in names}
		for cl in inClust:
			clustMemList.append(len(cl))
			myRepDictNo = {nL: 0 for nL in names}
			myRepDictLen = {nL: 0 for nL in names}
			maxNo = 0
			maxNoName = ''
			maxLen = 0
			maxLenName = ''
			totalLen = 0
			
			#get representations for each name
			for c in cl:
				mylen = int(c.rsplit('_',2)[1])
				for nL in names:
					#print c, nL, string.find(c, nL)
					if string.find(c,nL) != -1:
						myRepDictNo[nL] += 1
						myRepDictLen[nL] += mylen
						totalLen += mylen			
						break
							
			if len(cl) >= sizeThreshold:
				#print cl
				for nL in myRepDictNo:
					#print myRepDictNo[nL], maxNo
					if myRepDictNo[nL] > maxNo:
						maxNo = myRepDictNo[nL]
						maxNoName = nL
				#if GrDictNo[maxNoName] < maxNo:
				#	GrDictNo[maxNoName] = maxNo
				for nL in myRepDictLen:
					if myRepDictLen[nL] > maxLen:
						maxLen = myRepDictLen[nL]
						maxLenName = nL
				#if GrDictLen[maxLenName] < maxLen:
				#	GrDictLen[maxLenName] = maxLen
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
		ZAllNo = 0
		TPAllNo = 0
		FPAllNo = 0
		ZAllLen = 0
		TPAllLen = 0
		FPAllLen = 0
		for na in names:
			#if ZDict[na] and TPDict[na]:	  # tacoa p 13 "The overall specificity is computed over those classes that have a defined specificity value"
			#print na
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
		#print Zk
		if TPAllNo:
			SnAllNo = float(TPAllNo)/float(ZAllNo)
			SpAllNo = float(TPAllNo)/float(TPAllNo+FPAllNo)
		if TPAllLen:
			SnAllLen = float(TPAllLen)/float(ZAllLen)
			SpAllLen = float(TPAllLen)/float(TPAllLen+FPAllLen)
		F1No = 2*(SnAllNo*SpAllNo)/(SnAllNo+SpAllNo)
		F1Len = 2*(SnAllLen*SpAllLen)/(SnAllLen+SpAllLen)
		
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
		
		#nmi = NMI(inClust, corrClustReal)
		my_label = generateLabeling(inClust)
		idSubset = []
		for cl in inClust:
			for c in cl:
				idSubset.append(c)
		tr_label = generateLabelingSubset(corrClustReal,idSubset)
		nmi = sklearn.metrics.normalized_mutual_info_score(tr_label,my_label)
		ami = sklearn.metrics.adjusted_mutual_info_score(tr_label,my_label)
		vscore = sklearn.metrics.v_measure_score(tr_label,my_label)
		#ami = AMI(inClust, corrClustReal)
		
		avgClustSize = float(sum(clustMemList))/float(len(clustMemList))
		minClustSize = min(clustMemList)
		maxClustSize = max(clustMemList)
		numberClusters = 0
		for c in clustMemList:
		#	print c,
			if c >= sizeThreshold:
				numberClusters += 1
		
		#pprint.pprint(ZDictLen)
		
		#"Source;Score;Cut;N;J;L;SizeThreshold;NumberClusters;AvgClustSize;MinClustSize;MaxClustSize;NMI;AMI;V-score;SnAllNo;SpAllNo;RepFracNo;SnAllLen;SpAllLen;repFracLen"
		#pref = "NIL"
		outF.write("{!s};{!s};{!s};{:01.2f};{:01.2f};{:01.2f};{!s};".format(mText,score,cText,n,j,l,pref))
		outF.write("{!s};{!s};{:03.2f};{!s};{!s};{:01.4f};".format(sizeThreshold,numberClusters,avgClustSize,minClustSize,maxClustSize,numberClusters/33.0))
		outF.write("{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}".format(nmi,ami,vscore,SnAllNo,SpAllNo,F1No,repFracNo,SnAllLen,SpAllLen,F1Len,repFracLen))
		#for na in sorted(names):
		#	outF.write(";{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f};{:01.4f}".\
		#		format(TPDictNo[na],FPDictNo[na],GrDictNo[na],ZDictNo[na],TPDictLen[na],FPDictLen[na],GrDictLen[na],ZDictLen[na]))
		outF.write("\n")
		ctr += 1
		
	outF.close()
	print "done"
	
	
def generateCorrectClusteringFromResults(resPickle, nameFile, outFile, threshold):
	import cPickle
	import string
	nf = open(nameFile,'r')
	names = [l.rstrip() for l in nf.readlines()]
	corrDict = {n:[] for n in names}
	clustering = cPickle.load(open(resPickle,'rb'))
	for cluster in clustering:
		leaves = cluster.getLeaves()
		for contig in leaves:
			size = int(contig.rsplit('_',2)[-2])
			if size >= threshold:
				source = ''
				for possible in names:
					if string.find(possible,contig) != -1:
						source = possible
						corrDict[source].append(contig)
						break
	outList = corrDict.values()
	cPickle.dump(outList,open(outFile,'wb'))
	

def processDistanceInfo(inClusters, inDistance, nameFile):
	import cPickle as pickle
	import horatioUtils as hutil
	#associate roots - and index numbers - with max representation
	nf = open(nameFile,'r')
	nameList = [n.rstrip() for n in nf.readlines()]
	nf.close()
	clusters = pickle.load(open(inClusters,'rb'))
	rootDict = {clusters[c].root:'' for c in clusters.keys()}
	#print rootDict
	for c in clusters:
		_, maxName = purityOfCluster(clusters[c].getLeaves(),nameList)
		rootDict[clusters[c].root] = maxName
	distF = open(inDistance,'r')
	rootList = distF.readline().rstrip().split(",")
	rightDists = []
	wrongDists = []
	distF.readline()
	line = distF.readline().rstrip()
	ind = 0
	clustering = {n:[] for n in nameList}
	while line:
		myRoot = rootList[ind]
		myMax = rootDict[myRoot]
		dists = line.split(", ")
		clustering[myMax].append(ind)
		for j in range(len(dists)):
			dInfo = dists[j]
			[d,i] = dInfo.split(":")
			theirRoot = rootList[int(i)]
			theirMax = rootDict[theirRoot]
			if theirMax == myMax:
				rightDists.append(float(d))
			else:
				wrongDists.append(float(d))
		ind += 1
		line = distF.readline().rstrip()
	allDist = hutil.makeDistanceMatrix(inDistance)
	clustList = clustering.values()
	return rightDists,wrongDists,clustList,rootList,allDist
	

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


def processOrderedDistMatrix(orderedDistMatrix):
	from collections import Counter
	odmFile = open(orderedDistMatrix,'r')
	db_names = odmFile.readline().rstrip().split(',')
	db_genera = [n.split('_')[0] for n in db_names]
	contig_names = odmFile.readline().rstrip().split(',')
	contig_genera = [n.split('_')[0] for n in contig_names]
	contig_genera_uniq = sorted(list(set(contig_genera)))
	TPDict = {n:0 for n in contig_genera_uniq}
	FPDict = {n:0 for n in contig_genera_uniq}
	FNDict = {n:0 for n in contig_genera_uniq}
	ZDict = Counter(contig_genera)
	ct = 0
	uniqueClusters = {}
	ln = odmFile.readline().rstrip()
	while ln:
		entry = ln.split(',')[0]
		index = int(entry.split(':')[1])
		match_genus = db_genera[index]
		if match_genus in uniqueClusters.keys():
			uniqueClusters[match_genus] += 1
		else:
			uniqueClusters[match_genus] = 1
		contig_genus = contig_genera[ct]
		if contig_genus == match_genus:
			TPDict[contig_genus] += 1
		else:
			if match_genus in contig_genera_uniq:
				FPDict[match_genus] += 1
		ct += 1
		ln = odmFile.readline().rstrip()
	TPAll = 0
	FPAll = 0
	ZAll = 0
	for na in contig_genera_uniq:
		TPAll += TPDict[na]
		FPAll += FPDict[na]
		ZAll += ZDict[na]
	SN = float(TPAll)/float(ZAll)
	SP = float(TPAll)/float(TPAll+FPAll)
	print "SN: {:01.2f}".format(SN)
	print "SP: {:01.2f}".format(SP)
	print "NumClusters: {!s}".format(len(uniqueClusters.keys()))
	
	
def getPhyloDistribution(inClust,phyloDist,outPrefix):
	import cPickle
	import string
	import numpy
	import horatioUtils as hutil
	phyloF = open(phyloDist,'r')
	phyloLines = phyloF.readlines()
	print len(phyloLines)
	phyloNames = [string.replace(n," ","_") for n in phyloLines[0].rstrip().split(",")]
	phyloDists = []
	for x in phyloLines[1:]:
		ln = x.rstrip().split(",")
		phyloDists.append([float(n) for n in ln])
	print len(phyloDists), len(phyloDists[0])
	clustering = cPickle.load(open(inClust,'rb'))
	ct = 0
	orderingOfSources = {}
	myIntervals = []
	for cluster in clustering:
		intervalStart = ct
		contigs = cluster.getLeaves()
		for con in contigs:
			for pn in phyloNames:
				if string.find(con,pn) != -1:
					orderingOfSources[ct] = pn
					break
			ct += 1
		intervalEnd = ct
		myIntervals.append(range(intervalStart,intervalEnd))
	numContigs = ct
	distances = -1.0 * numpy.ones((numContigs,numContigs))
	for x in range(numContigs):
		distances[x][x] = 0.0
		srcX = orderingOfSources[x]
		indX = phyloNames.index(srcX)
		for y in range(x):
			srcY = orderingOfSources[y]
			indY = phyloNames.index(srcY)
			dist = phyloDists[indX][indY]
			distances[x][y] = dist
			distances[y][x] = dist
	hutil.csv2DistanceDistribution(distances,outPrefix+".intra.csv",numContigs,intervals=myIntervals,include=1)
	hutil.csv2DistanceDistribution(distances,outPrefix+".inter.csv",numContigs,intervals=myIntervals,include=0)

