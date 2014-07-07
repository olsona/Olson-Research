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
	print graph
	print start
	if visited is None:
		visited = set()
	print visited
	visited.add(start)
	for next in graph[start] - visited:
		dfs(graph, next, visited)
	return visited			
	
	
def makeCorrectClustering(contigFile, nameFile, out, threshold=0):
	import string
	import cPickle as pickle
	cf = open(contigFile,'r')
	nf = open(nameFile,'r')
	nameList = []
	for l in nf.readlines():
		n = l.rstrip()
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