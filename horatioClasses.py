import horatioUtils as hutil
#from bootstrapCorrectness import *

class Cluster:
    def __init__(self, seed, dict = None, root = None, closeList = None):
	self.seed = seed
	if root is None:
	    self.root = self.seed
	else:
	    self.root = root	
	if dict is None:
	    self.dict = {self.root:[]}
	else:
	    self.dict = dict	
	if closeList is None:
	    self.closeList = {}
	else:
	    self.closeList = closeList
			
    def __str__(self):
	return "Seed: {!s}\tLeaves: {!s}".format(self.seed,self.getLeaves())

    def getLeaves(self):
	if self.root is None or self.dict is None:
	    return []
	else:
	    return hutil.getLeavesUtil(self.dict, self.root)
    
    def addNode(self, parent, child):
	if parent in self.dict:
	    self.dict[parent].append(child)
	else:
	    self.dict[parent] = [child]


    def addMatch(self, goodMatch):
	name = goodMatch[0]
	score = goodMatch[1]
	if name in self.closeList:
	    self.closeList[name].append(score)
	else:
	    self.closeList[name] = [score]
			
    def getMostCommonNeighbor(self):
	total = 0
	max = 0
	maxName = ''
	for n in self.closeList:
	    l = len(self.closeList[n])
	    total += l
	    if l > max:
	        max = l
	        maxName = n
	return float(max)/float(total), maxName
		
    def getNeighborInfo(self):
	total = 0
	repDict = {n:0 for n in self.closeList}
	for n in self.closeList:
	    l = len(self.closeList[n])
	    total += l
	    repDict[n] = l
	for n in repDict:
	    tmp = repDict[n]
	    repDict[n] = float(tmp)/float(total)
	    sortTuple = sorted(repDict.items(), key=lambda x: x[1])
	return sortTuple[::-1]

    def addClusters(self, otherClusters, ubercontig):
	subs = [self.root]
	for o in otherClusters:
	    subs.append(o.root)
	self.dict.update(o.dict)
	self.root = ubercontig
	self.dict[ubercontig] = subs
	
    # ***
    #def purityMax(self, names):
    #	if self.dict is None:
    #	    return None, None
    #	else:
    #	    return purityOfCluster(self.getAll(), names)
    # ***
		
# THE TREES ARE UPSIDE DOWN!  dict[j] is all of the *parents* of j!


class Contig:
    def __init__(self, name, goodMatches=None):
	self.name = name
	if goodMatches is None:
	    self.goodMatches = []
	else:
	    self.goodMatches = goodMatches

    def __str__(self):
	return "Name: {!s}".format(self.name)