from bootstrapUtils import *

class Cluster:
	def __init__(self, seed, dict = None, root = None, closeList = None):
		self.seed = seed
		if dict is None:
			self.dict = {}
		else:
			self.dict = dict
		if root is None:
			self.root = self.seed
		else:
			self.root = root
		if closeList is None:
			self.closeList = {}
		else:
			self.closeList = closeList
		#print self.seed

	def __str__(self):
		return "Seed: {!s}\tLeaves: {!s}".format(self.seed,self.get_leaves())

	def getLeaves(self):
		import pprint
		pprint.pprint(self.dict)
		if self.root is None:
			return []
		elif not any(self.dict):
			return []
		else:
			return getLeaves(self.dict, self.root)

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
			return self.getLeaves().append(self.root)

	def addMatch(self, goodMatch):
		name = goodMatch[0]
		score = goodMatch[1]
		if name in self.closeList:
			self.closeList[name].append(score)
		else:
			self.closeList[name] = [score]
			
	def getMatches(self, nameList):
		import string
		matchDict = {name:0 for name in nameList}
		total = 0
		for gm in self.closeList:
			for n in nameList:
				if string.find(gm, n) != -1:
					matchDict[n] += 1
					total += 1
		#http://stackoverflow.com/a/2258273
		sorted = sorted(matchDict.items(), key=lambda x: x[1])
		return sorted[::-1], total

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