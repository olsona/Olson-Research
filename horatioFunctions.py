import os
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
	
	
def scoreRAIphyFinal(DB, fSeed, toMatch, computePath, outputFile):
    os.system("{!s}rait -new -i {!s} -o {!s} >/dev/null 2>&1".format(computePath, fSeed, DB))
    os.system("{!s}rai -I {!s} -d {!s} >/dev/null 2>&1".format(computePath, toMatch, DB))
    short = toMatch.rsplit("/",1)[1]
    os.system("cp {!s}/{!s}.bin {!s}".format(os.getcwd(), short, outputFile+"_dists_sorted")) # moves results to results folder
    os.system("rm {!s}/{!s}.bin".format(os.getcwd(), short))
	
	
def scoreTETRA(DB, fSeed, matches, toMatch, allContigs):
    f = open(toMatch+"-2",'r')
    for l in f.readlines():
        sp = l.rstrip().split("\t")
	nm = sp[0]
	co = Contig(nm)
	allContigs[nm] = co
    mDB = "{!s}_M".format(DB)
    os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s} >/dev/null".format(fSeed,DB))
    os.system("perl tetraZscores.pl -k 4 -m {!s}-2 {!s} >/dev/null".format(toMatch,mDB))
    os.system("perl tetraCorrelation.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,matches))
    #os.system("rm {!s}".format(mDB))
	
	
def scoreTETRAFinal(DB, fSeed, toMatch, outputFile):
    os.system("perl tetraZscores.pl -k 4 -m {!s} {!s} >/dev/null".format(fSeed,DB))
    os.system("perl tetraCorrelation.pl {!s} {!s} {!s}_dists_sorted >/dev/null".format(DB, DB, outputFile))
	
	
def scoreTACOA(DB, fSeed, matches, toMatch, allContigs):
    f = open(toMatch+"-2",'r')
    for l in f.readlines():
	sp = l.rstrip().split("\t")
	nm = sp[0]
	co = Contig(nm)
	allContigs[nm] = co
    mDB = "{!s}_M".format(DB)

    #print("  Running1: perl tacoaCount.pl -k 4 {!s}-2 {!s}".format(fSeed,DB))
    #d1 = os.path.getsize(fSeed+"-2")
    #os.system("du -h {!s}-2".format(fSeed))
    os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(fSeed,DB))
    
    #print("  Running2: perl tacoaCount.pl -k 4 {!s}-2 {!s}".format(toMatch,mDB))
    #d2 = os.path.getsize(toMatch+"-2")
    #os.system("du -h {!s}-2".format(toMatch))
    os.system("perl tacoaCount.pl -k 4 {!s}-2 {!s} >/dev/null".format(toMatch,mDB))

    #print("  Running3: perl tacoaDistance.pl {!s} {!s} {!s}".format(DB,mDB,matches))
    #d3 = os.path.getsize(DB)
    #d4 = os.path.getsize(mDB)
    #os.system("du -h {!s} {!s}".format(DB, mDB))
    os.system("perl tacoaDistance.pl {!s} {!s} {!s} >/dev/null".format(DB,mDB,matches))
    #print "{:2.2e} {:2.2e} {:2.2e} {:2.2e}".format(d1,d2,d3,d4)
	
	
def scoreTACOAFinal(DB, fSeed, toMatch, outputFile):
    os.system("perl tacoaCount.pl -k 4 {!s} {!s} >/dev/null".format(fSeed,DB))
    os.system("perl tacoaDistance.pl {!s} {!s} {!s}_dists_sorted >/dev/null".format(DB, DB, outputFile))