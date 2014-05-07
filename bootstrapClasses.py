from bootstrapUtils import *
from bootstrapCorrectness import *

class Cluster:
	def __init__(self, seed, dict = None, root = None, closeListSeed = None, closeListMax = None):
		self.seed = seed
		
		if root is None:
			self.root = self.seed
		else:
			self.root = root
			
		if dict is None:
			self.dict = {self.root:[]}
		else:
			self.dict = dict
			
		if closeListSeed is None:
			self.closeListSeed = {}
		else:
			self.closeListSeed = closeListSeed
			
		if closeListMax is None:
			self.closeListMax = {}
		else:
			self.closeListMax = closeListMax

	def __str__(self):
		return "Seed: {!s}\tLeaves: {!s}".format(self.seed,self.getLeaves())

	def getLeaves(self):
		import pprint
		if self.root is None:
			return []
		else:
			return getLeavesUtil(self.dict, self.root)

	def purityMax(self, names):
		if self.dict is None:
			return None, None
		else:
			return purityOfCluster(self.getAll(), names)

	def addNode(self, parent, child):
		if parent in self.dict:
			self.dict[parent].append(child)
		else:
			self.dict[parent] = [child]

	def getAll(self):
		if self.root is None or self.dict is None:
			return []
		else:
			res = self.getLeaves()
			return res

	def addMatchSeed(self, goodMatch):
		name = goodMatch[0]
		score = goodMatch[1]
		if name in self.closeListSeed:
			self.closeListSeed[name].append(score)
		else:
			self.closeListSeed[name] = [score]
			
	def addMatchMax(self, goodMatch):
		name = goodMatch[0]
		score = goodMatch[1]
		if name in self.closeListMax:
			self.closeListMax[name].append(score)
		else:
			self.closeListMax[name] = [score]
			
	def getMostCommonNeighbor(self):
		total = 0
		max = 0
		maxName = ''
		for n in self.closeListSeed:
			l = len(self.closeListSeed[n])
			total += l
			if l > max:
				max = l
				maxName = n
		return float(max)/float(total), maxName
		
	def getNeighborInfo(self):
		total = 0
		repDict = {n:0 for n in self.closeListSeed}
		for n in self.closeListSeed:
			l = len(self.closeListSeed[n])
			total += l
			repDict[n] = l
		for n in repDict:
			tmp = repDict[n]
			repDict[n] = float(tmp)/float(total)
		sortTuple = sorted(repDict.items(), key=lambda x: x[1])
		return sortTuple[::-1]

	def addCluster(self, others, ubercontig):
		subs = [self.root]
		for o in others:
			subs.append(o.root)
			self.dict.update(o.dict)
		self.root = ubercontig
		self.dict[ubercontig] = subs
		
	
		

# THE TREES ARE UPSIDE DOWN!  dict[j] is all of the *parents* of j!


class Contig:
	def __init__(self, name, myCluster=None, goodMatchesSeed=None, goodMatchesMax=None):
		self.name = name
		if myCluster is None:
			self.myCluster = None
		else:
			self.myCluster = myCluster
		if goodMatchesSeed is None:
			self.goodMatchesSeed = []
		else:
			self.goodMatchesSeed = goodMatchesSeed
		if goodMatchesMax is None:
			self.goodMatchesMax = []
		else:
			self.goodMatchesMax = goodMatchesMax

	def __str__(self):
		return "Name: {!s}\nFile: {!s}\nCluster: {!s}".format(self.name, self.myCluster)