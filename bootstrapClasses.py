from bootstrapUtils import *

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
		if self.root is None:
			return []
		else:
			return getLeavesUtil(self.dict, self.root)

	def purityMax(self, names):
		if self.dict is None:
			return None, None
		else:
			return purityOfCluster(dict, names)

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

	def addMatch(self, goodMatch):
		name = goodMatch[0]
		score = goodMatch[1]
		if name in self.closeList:
			self.closeList[name].append(score)
		else:
			self.closeList[name] = [score]

	def addCluster(self, others, ubercontig):
		subs = [self.root]
		for o in others:
			subs.append(o.root)
			self.dict.update(o.dict)
		self.root = ubercontig
		self.dict[ubercontig] = subs
		
	
		

# THE TREES ARE UPSIDE DOWN!  dict[j] is all of the *parents* of j!


class Contig:
	def __init__(self, name, myCluster=None, goodMatches=None):
		self.name = name
		if myCluster is None:
			self.myCluster = None
		else:
			self.myCluster = myCluster
		if goodMatches is None:
			self.goodMatches = []
		else:
			self.goodMatches = goodMatches

	def __str__(self):
		return "Name: {!s}\nFile: {!s}\nCluster: {!s}".format(self.name, self.myCluster)