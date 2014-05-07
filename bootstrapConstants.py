from bootstrapUtils import *
from bootstrapFunctions import *
from bootstrapCorrectness import *

myMethods = ['raiphy']

usageString = 'bootstrap.py -i <input file> -o <output file> -p <score computation path> '\
			'[-c <cut schedule> -s <scoring function>]'
				
bigUsageString = '-i, --ifile\t\tInput metagenomics FASTA file\n'+\
				'-o, --ofile\t\tOutput file\n'+\
				'-p, --path\t\tPath to extra code for computation\n'+\
				'-s, --score\t\tScoring function for use (optional)\n'+\
				'-c, --cut\t\tPartition schedule (optional)\n'+\
				'\t\t\t	 format: [1,2,5]\n'+\
				'\t\t\t	 [1,2,5] implies the partition schedule > 5kbp, 2kbp, 1kbp\n'
				
defaultSchedule = [2,6,10,14,18,22]

#http://my.safaribooksonline.com/book/programming/python/0596007973/python-shortcuts/pythoncook2-chp-4-sect-16
correctnessDictSeed = {
	"species":	checkCorrectSpeciesOlsonFormat,
	"genus":	checkCorrectGenusOlsonFormat,
}

correctnessDictMax = {
	"species":	checkCorrectMatchClusterMaxSpecies,
	"genus":	checkCorrectMatchClusterMaxGenus,
}

scoringMethod = {
	"raiphy":	scoreRAIphy,
	"tetra":	scoreTETRA,
	"tacoa":	scoreTACOA,
}

scoringMethodFinal = {
	"raiphy":	scoreRAIphyFinal,
	"tetra":	scoreTETRAFinal,
	"tacoa":	scoreTACOAFinal,
}