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
	f.close()
	return seq_name, concat


def makeDistanceMatrix(scoreFile):
	f = open(scoreFile,'r')
	lns = f.readlines()
	rows = len(lns)-2
	cols = len(lns[2].split(", "))
	dists = [[0.0 for c in range(cols)] for r in range(rows)]
	for r in range(rows):
		li = lns[r+2].split(", ")
		for l in li:
			[sStr, cStr] = l.split(":")
			s = float(sStr)
			c = int(cStr)
			dists[r][c] = s
	return dists


def getMatchesAboveThreshold(inline, threshold):
	ln = inline.rstrip().split(", ")
	out = []
	for i in range(len(ln)):
		info = ln[i].split(":")
		match = int(info[0])
		dist = float(info[1])
		if dist >= threshold:
			out.append([match,dist])
		else:
			break
	return out


def getMatchesWithinPercentage(inline, pct):
	ln = inline.rstrip().split(", ")
	out = []
	#find best score
	info = ln[0].split(":")
	match = int(info[0])
	best = float(info[1])
	out.append([match,best])
	threshold = best*(1.0-float(pct)/100.0)
	#find all scores within pct of
	for i in range(1,len(ln)):
		info = ln[i].split(":")
		match = int(info[0])
		dist = float(info[1])
		if dist >= threshold:
			out.append([match,dist])
		else:
			break
	return out
	

def dfs(graph, start, visited=None):
	# http://eddmann.com/posts/depth-first-search-and-breadth-first-search-in-python/
	if visited is None:
		visited = set()
	visited.add(start)
	for next in graph[start] - visited:
		dfs(graph, next, visited)
	return visited			
	
	
def processAPLabels(labels, items):
	lDict = {}
	for i in range(len(labels)):
		l = labels[i]
		if l in lDict.keys():
			lDict[l].append(items[i])
		else:
			lDict[l] = [items[i]]
	return lDict.values()
	
	
def makeCorrectClustering(contigFile, nameFile, out, threshold=0):
	import string
	import cPickle as pickle
	cf = open(contigFile,'r')
	nf = open(nameFile,'r')
	nameList = []
	for l in nf.readlines():
		n = l.rstrip().split('_')[0]
		if len(n) > 0:
			nameList.append(n)
	resDict = {name:[] for name in nameList}
	l = cf.readline()
	while l:
		if l[0] == ">":
			#we have a contig name
			li = l.rstrip()[1:]
			src = ''
			for n in nameList:
				if string.find(li, n) != -1:
					src = n
					break
			leng = int(li.split("_")[-2])
			if leng >= threshold:
				 resDict[src].append(li)
		l = cf.readline()
	nf.close()
	cf.close()
	resList = []
	for r in resDict:
		resList.append(resDict[r])
	pickle.dump(resList, open(out,"wb"))
	

def percentile90(array):
	import numpy
	return numpy.percentile(array,90)


def percentile80(array):
	import numpy
	return numpy.percentile(array,80)	


def percentile70(array):
	import numpy
	return numpy.percentile(array,70)
	

def percentile60(array):
	import numpy
	return numpy.percentile(array,60)
	

def percentile40(array):
	import numpy
	return numpy.percentile(array,40)

	
def makeGenomeDistanceMatrix(namesFile, names2IDFile, distFile, outCSV):
	import numpy
	import subprocess
	import string
	import itertools
	# populate names file
	nF = open(namesFile,'r')
	names = [string.replace(n.rstrip(),'_',' ') for n in nF.readlines()]
	N = len(names)
	distMatrix = numpy.zeros((N,N)) - 1.0
	names2IDsDict = {n:[] for n in names}
	# map names to IDs
	for n in names:
		[genus,species,_] = n.split(" ",2)
		result = subprocess.check_output("grep {ge} {fi} | grep {sp}".format(ge=genus,fi=names2IDFile,sp=species),shell=True).split("\n")
		for r in result:
			id = r.split("\t")[0]
			if len(id) > 0:
				names2IDsDict[n].append(id)
	# populate distance matrix
	for x in range(len(names)):
		n1 = names[x]
		index1 = names.index(n1)
		try:
			idSet1 = names2IDsDict[n1]
			idRegex1 = "\|".join(idSet1)
			subprocess.check_output("grep '{i1}' {fi} > {fi}.current".format(i1=idRegex1,fi=distFile),shell=True)
			distMatrix[x][x] = 0.0
			for y in range(x):
				print x, y
				n2 = names[y]
				index2 = names.index(n2)
				if n1 != n2:
					# get distance between distinct species
					idSet2 = names2IDsDict[n2]
					idRegex2 = "\|".join(idSet2)
					try:
						result = subprocess.check_output("grep '{i2}' {fi}.current".format(i2=idRegex2,fi=distFile),shell=True)
						dist = float(result.split("\n")[0].split("\t")[2])
						distMatrix[index1][index2] = dist
						distMatrix[index2][index1] = dist
						print "{0}, {1}, {2}".format(n1,n2,dist)
					except subprocess.CalledProcessError:
						print "Error"
						pass
				else:
					distMatrix[index1][index2] = 0.0
		except subprocess.CalledProcessError:
			distMatrix[index1][index1] = 0.0
	return names, distMatrix
	

def csv2DistanceDistribution(distMatrix, out, num, intervals = [], include = 1):
    # http://stackoverflow.com/a/5328669
    import numpy as np
    
    mylist = []
    if intervals == []:
        for i in range(len(distMatrix)):
            row = distMatrix[i]
            for j in range(i):
                if row[j] > 0:
                    mylist.append(row[j])
    else:
        if include:
            for n in range(len(intervals)):
                for i in intervals[n]:
                    row = distMatrix[i]
                    for j in set(intervals[n]) - set([i]):
                        mylist.append(row[j])
        else:
            for n in range(len(intervals)):
                for i in intervals[n]:
                    row = distMatrix[i]
                    okset = set(range(num)) - set(intervals[n])
                    for j in okset:
                        mylist.append(row[j])

    mat = np.array(mylist)
    np.savetxt(out, mat, fmt='%.5e', delimiter=',')