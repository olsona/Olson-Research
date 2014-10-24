def removeSomeKnownSpecies(dbFolder, keyNames, percentKeep, outList):
	import glob
	import random
	backgroundList = []
	keyList = []
	listFiles = glob.glob(dbFolder + "/*.fna")
	for fi in listFiles:
		if isInName(fi, keyNames):
			keyList.append(fi)
		else:
			backgroundList.append(fi)
	nKeys = len(keyList)
	nKeep = int((percentKeep * nKeys)+0.5)
	keepList = random.sample(keyList, nKeep)
	outF = open(outList, 'w')
	for bL in backgroundList:
		name = bL.split("/")[-1].split(".")[0]
		outF.write(name + "\t" + bL + "\n")
	for kL in keepList:
		name = kL.split("/")[-1].split(".")[0]
		outF.write(name + "\t" + kL + "\n")
	outF.close()
		

def isInNameSet(fileName, nameList):
	import string
	res = False
	for n in nameList:
		if string.find(fileName, n) > -1:
			res = True
			break
	return res
	
	
def doTest(score, path, metagenome, dbFolder, names, percentKeep, workFolder, out):
	import os
	import horatioFunctions as hfun
	removeSomeKnownSpecies(dbFolder, names, percentKeep, workFolder + "/DB")
	os.system("mkdir {!s}/contigs/".format(workFolder))
	os.system("perl sepMetagenome.pl {!s}/contigs/ {!s} mtgnm".format(workFolder, metagenome))
	if score == "raiphy":
		hfun.scoreRAIphyFinal(workFolder+"/DB", workFolder + "mtgnm", path, out)
	elif score == "tetra":
		hfun.scoreTETRAFinal(workFolder+"/DB", workFolder + "mtgnm", out)
	elif score == "tacoa":
		hfun.scoreTACOAFinal(workFolder+"/DB", workFolder + "mtgnm", out)