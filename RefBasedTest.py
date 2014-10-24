import os

def removeSomeKnownSpecies(dbFolder, keyNames, percentKeep, outList):
	import glob
	import random
	backgroundList = []
	keyList = []
	listFiles = glob.glob(dbFolder + "/*.fna")
	for fi in listFiles:
		if isInNameSet(fi, keyNames):
			keyList.append(fi)
		else:
			backgroundList.append(fi)
	nKeys = len(keyList)
	nKeep = int(((percentKeep/100.0) * nKeys)+0.5)
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


def makeRaiDB(path, list, outDB):
	os.system("{!s}rait -new -i {!s} -o {!s} >/dev/null 2>&1".format(path, list, outDB))
	
def makeTetraDB(list, outDB):
	os.system("perl tetraZscores.pl -k 4 -m {!s} {!s} >/dev/null".format(list,outDB))
	
def makeTacoaDB(list, outDB):
	os.system("perl tacoaCount.pl -k 4 {!s} {!s} >/dev/null".format(list,outDB))
	
def scoreRai(path, list, DB, outputFile):
	os.system("{!s}rai -I {!s}-1 -d {!s} >/dev/null 2>&1".format(path, list, DB))
	short = toMatch.rsplit("/",1)[1]
	os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, outputFile)) # moves results to results folder
	os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))

def scoreTetra(list, DB, outputFile):
	mDB = DB + ".m"
	os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s}".format(list, mDB))
	os.system("perl tetraCorrelation.pl {!s} {!s} {!s} >/dev/null".format(DB, mDB, outputFile))
	
def scoreTacoa(list, DB, outputFile):
	mDB = DB + ".m"
	os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(list,mDB))
	os.system("perl tacoaDistance.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,outputFile))
	
	
def doTest(score, path, metagenomes, suffixes, dbFolder, names, percentKeep, workFolder, numDB, out):
	import os
	import horatioFunctions as hfun
	os.system("mkdir {!s}/contigs/".format(workFolder))
	for i in range(numDB):
		myDList = workFolder + "/list." + str(i)
		removeSomeKnownSpecies(dbFolder, names, percentKeep, myDList)
		myDB = workFolder + "/DB." + str(i)
		if score == "raiphy":
			makeRaiDB(path, myDList, myDB)
		elif score == "tacoa":
			makeTacoaDB(myDList, myDB)
		elif score == "tetra":
			makeTetraDB(myDList, myDB)
		print "Made DB"
		for j in range(len(metagenomes)):
			mg = metagenomes[j]
			os.system("perl fasta2tab.pl {!s} {!s}/mtgnm.tab".format(mg, workFolder))
			os.system("perl sepMetagenome.pl {!s}/contigs/ {!s}/mtgnm.tab {!s}/mtgnm".format(workFolder, workFolder, workFolder))
			list = workFolder + "mtgnm"
			myOut = out
			if score == "raiphy":
				scoreRai(path, list, myDB, out + "." + suffixes[j])
			elif score == "tetra":
				scoreTetra(list, myDB, out + "." + suffixes[j])
			elif score == "tacoa":
				scoreTacoa(list, myDB, out + "." + suffixes[j])
			print i, j
			
			
	