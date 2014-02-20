from bootstrapUtils import *

class Cluster:
    def __init__(self, seed, dict = None, root = None):
        self.seed = seed
        if dict is None:
            self.dict = {}
        else:
            self.dict = dict
        if root is None:
            self.root = None
        else:
            self.root = root

    def __str__(self):
        return "Seed: {!s}\nLeaves: {!s}".format(self.seed,self.get_leaves())

    def get_leaves(self):
        if self.root is None or self.dict is None:
            return []
        else:
            return getLeaves(self.dict, self.root)

    def purityMax(self, names):
        if self.dict is None:
            return None, None
        else:
            return purityOfCluster(dict, names)

    def addNode(self, parent, child):
        print "I am trying to add a node"
        if parent in self.dict:
            self.dict[parent].append(child)
        else:
            self.dict[parent] = [child]


# THE TREES ARE UPSIDE DOWN!  dict[j] is all of the *parents* of j!


class Contig:
    def __init__(self, name, file, myCluster=None, good_matches=None):
        self.name = name
        self.file = file
        if myCluster is None:
            self.myCluster = None
        else:
            self.myCluster = myCluster
        if good_matches is None:
            self.good_matches = []
        else:
            self.good_matches = good_matches

    def __str__(self):
        return "Name: {!s}\nFile: {!s}\nCluster: {!s}".format(self.name, self.file, self.myCluster)