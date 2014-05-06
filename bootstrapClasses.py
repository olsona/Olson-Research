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
			print "Root:", self.root
			print "Dict:"
			pprint.pprint(self.dict)
			res = getLeavesUtil(self.dict, self.root)
			print res
			print "Done\n"
			return res

	def purityMax(self, names):
		if self.dict is None:
			return None, None
		else:
			print self.getLeaves(), self.seed
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
			res.append(self.root)
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