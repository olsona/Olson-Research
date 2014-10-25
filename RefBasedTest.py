import os

def removeSomeKnownSpecies(dbFolder, keyNames, percentKeep, outList):
	import glob
	import random
	keyList = []
	listFiles = glob.glob(dbFolder + "/*.fna")
	for fi in listFiles:
		if isInNameSet(fi, keyNames):
			keyList.append(fi)
	nKeys = len(keyList)
	nKeep = int(((percentKeep/100.0) * nKeys)+0.5)
	keepList = random.sample(keyList, nKeep)
	outF = open(outList, 'w')
	for kL in keepList:
		name = kL.split("/")[-1].split(".")[0]
		outF.write(name + "\t" + kL + "\n")
	outF.close()
		

def makeBackgroundList(dbFolder, keyNames, outList):
	import glob
	backgroundList = []
	listFiles = glob.glob(dbFolder + "/*.fna")
	for fi in listFiles:
		if isInNameSet(fi, keyNames):
			pass
		else:
			backgroundList.append(fi)
	outF = open(outList, 'w')
	for bL in backgroundList:
		name = bL.split("/")[-1].split(".")[0]
		outF.write(name + "\t" + bL + "\n")
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
	print list, mDB
	os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s}".format(list, mDB))
	print "made mDB"
	os.system("perl tetraCorrelation.pl {!s} {!s} {!s} >/dev/null".format(DB, mDB, outputFile))
	
def scoreTacoa(list, DB, outputFile):
	mDB = DB + ".m"
	os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(list,mDB))
	os.system("perl tacoaDistance.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,outputFile))
	
	
def doTest(path, metagenomes, suffixes, dbFolder, names, percentKeep, workFolder, numDB, out):
	import os
	import horatioFunctions as hfun
	os.system("mkdir {!s}/contigs/".format(workFolder))
	with open(names) as nF:
		namesList = [l.rstrip() for l in list(nF)]
	baseList = workFolder + "/base.list"
	makeBackgroundList(dbFolder, namesList, baseList)
	baseDB = workFolder + "/Base.db"
	makeRaiDB(path, baseList, baseDB + ".raiphy")
	print "Made base Raiphy DB"
	makeTacoaDB(baseList, baseDB + ".tacoa")
	print "Made base Tacoa DB"
	makeTetraDB(baseList, baseDB + ".tetra")
	print "Made base Tetra DB"
	for i in range(numDB):
		myDList = workFolder + "/list." + str(i)
		removeSomeKnownSpecies(dbFolder, namesList, percentKeep, myDList)
		myDB = workFolder + "/DB." + str(i)
		makeRaiDB(path, myDList, myDB+".ra")
		makeTacoaDB(myDList, myDB+".ta")
		makeTetraDB(myDList, myDB+".te")
		os.system("tail -n +1 {!s} > {!s}".format(myDB+".ra",myDB+".rai"))
		os.system("cat {!s}.raiphy {!s}.rai > {!s}.raiphy".format(baseDB, myDB, myDB))
		os.system("cat {!s}.tacoa {!s}.ta > {!s}.tacoa".format(baseDB, myDB, myDB))
		os.system("cat {!s}.tetra {!s}.te > {!s}.tetra".format(baseDB, myDB, myDB))
		print "Made DBs - Step " + str(i)
		for j in range(len(metagenomes)):
			mg = metagenomes[j]
			os.system("perl fasta2tab.pl {!s} {!s}/mtgnm.tab".format(mg, workFolder))
			os.system("perl sepMetagenome.pl {!s}/contigs/ {!s}/mtgnm.tab {!s}/mtgnm".format(workFolder, workFolder, workFolder))
			mtgList = workFolder + "/mtgnm"
			myOut = out + "." + str(i) + "." + suffixes[j]
			scoreRai(path, mtgList, myDB+".raiphy", myOut + ".raiphy")
			scoreTetra(mtgList, myDB+".tetra", myOut + ".tetra")
			scoreTacoa(mtgList, myDB+".tacoa", myOut + ".tacoa")
			print i, j
			
			
	