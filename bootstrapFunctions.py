import os
from bootstrapClasses import Contig

def scoreRAIphy(DB, raiPath, fSeed, matches, toMatch, allContigs):
	# os.system("{!s}rait -new -i {!s}-2 -o {!s} >/dev/null 2>&1".format(raiPath, fSeed, DB))
	os.system("{!s}rait -i {!s}-2 -o {!s} >/dev/null 2>&1".format(raiPath, fSeed, DB))
	# Process contigs to match
	f = open(toMatch+"-2",'r')
	for l in f.readlines():
		sp = l.rstrip().split("\t")
		nm = sp[0]
		fi = sp[1]
		co = Contig(nm,fi)
		allContigs[nm] = co
	# Match ith contigs to DB
	os.system("{!s}rai -I {!s}-1 -d {!s} >/dev/null 2>&1".format(raiPath, toMatch, DB))
	short = toMatch.rsplit("/",1)[1]
	os.system("cp {!s}/{!s}-1.bin {!s}".format(os.getcwd(), short, matches)) # moves results to results folder
	os.system("rm {!s}/{!s}-1.bin".format(os.getcwd(), short))
	
def scoreTETRA(fSeed, DB, matches, mDB):
    os.system("perl count.pl -k 4 -mf {!s}-2 {!s}".format(fSeed,DB))
    os.system("perl count.pl -k 4 -mf {!s} {!s}".format(matches,mDB))
    