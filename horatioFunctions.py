import os, sys
from horatioClasses import Contig

def scoreRAIphy(DB, raiPath, fSeed, matches, toMatch, allContigs):
	os.system("{!s}rait -i {!s}-2 -o {!s} >/dev/null 2>&1".format(raiPath, fSeed, DB))
	#os.system("{!s}rait -i {!s}-2 -o {!s}".format(raiPath, fSeed, DB))
	# Process contigs to match
	f = open(toMatch+"-2",'r')
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]
		co = Contig(nm)
		allContigs[nm] = co
	# Match ith contigs to DB
	os.system("{!s}rai -I {!s}-1 -d {!s} >/dev/null 2>&1".format(raiPath, toMatch, DB))
	short = toMatch.rsplit("/",1)[1]
	os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, matches)) # moves results to results folder
	os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))
	
	
def scoreRAIphyFinal(DB, fSeed, computePath, outputFile):
	os.system("{!s}rait -new -i {!s}-2 -o {!s} >/dev/null 2>&1".format(computePath, fSeed, DB))
	os.system("{!s}rai -I {!s}-1 -d {!s} >/dev/null 2>&1".format(computePath, fSeed, DB))
	short = fSeed.rsplit("/",1)[1]
	os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, outputFile+"_dists_sorted")) # moves results to results folder
	os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))
	
	
def scoreTETRA(DB, fSeed, matches, toMatch, allContigs):
	f = open(toMatch+"-2",'r')
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]
		co = Contig(nm)
		allContigs[nm] = co
	mDB = "{!s}_M".format(DB)
	os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s}".format(fSeed,DB)) 
	os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s}".format(toMatch,mDB))
	os.system("perl tetraCorrelation.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,matches))
	#os.system("rm {!s}".format(mDB))
	
	
def scoreTETRAFinal(DB, fSeed, outputFile):
	os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s} >/dev/null".format(fSeed,DB))
	os.system("perl tetraCorrelation.pl {!s} {!s} {!s}_dists_sorted >/dev/null".format(DB, DB, outputFile))
	
	
def scoreTACOA(DB, fSeed, matches, toMatch, allContigs):
	f = open(toMatch+"-2",'r')
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]
		co = Contig(nm)
		allContigs[nm] = co
	mDB = "{!s}_M".format(DB)
	os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(fSeed,DB))
	os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(toMatch,mDB))
	os.system("perl tacoaDistance.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,matches))
	
	
def scoreTACOAFinal(DB, fSeed, outputFile):
	#os.system("cat {!s}-2 | wc -l".format(fSeed))
	#os.system("head {!s}-2".format(fSeed))
	os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(fSeed,DB))
	#os.system("cat {!s} | wc -l".format(DB))
	#os.system("head {!s}".format(DB))
	os.system("perl tacoaDistance.pl {!s} {!s} {!s}_dists_sorted >/dev/null".format(DB, DB, outputFile))