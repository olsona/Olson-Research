from bootstrapUtils import *

class Cluster:
    def __init__(self, root, tree=None):
        self.root = root
        if tree is None:
            self.tree = {}
            self.tree[self.root] = []
        else:
            self.tree = dict

    def __str__(self):
        leaves = [self.root]
        leaves.append(self.get_leaves())
        return "Root: {!s}\nCluster: {!s}".format(self.root,leaves)

    def get_leaves(self):
        return getLeaves(self.tree, self.root)

    def purityMax(self, names):
        return purityOfCluster(tree, names)

    def addNode(self, parent, child):
        if parent in self.tree:
            self.tree[parent].append(child)
        else:
            self.tree[parent] = [child]


class Contig:
    def __init__(self, name, file, cluster=None, good_matches=None):
        self.name = name
        self.file = file
        if cluster is None:
            self.cluster = None
        else:
            self.cluster = cluster
        if good_matches is None:
            self.good_matches = []
        else:
            self.good_matches = good_matches