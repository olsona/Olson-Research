from bootstrapUtils import *
from bootstrapFunctions import *
from bootstrapCorrectness import *

myMethods = ['raiphy']

usageString = 'bootstrap.py -i <input file> -o <output file> -p <RAIphy path> '\
			'[-r <reference database> -c <cooling schedule>]'
				
bigUsageString = '-i, --ifile\t\tInput metagenomics FASTA file\n'+\
				'-o, --ofile\t\tOutput file\n'+\
				'-p, --partition\tpartition schedule (optional)\n'+\
				'\t\t\t	 format: [1,2,5]\n'+\
				'\t\t\t	 [1,2,5] implies the partition schedule > 5kbp, 2kbp, 1kbp\n'+\
				'-r, --raipath\t\tPath to RAIphy'
				
defaultSchedule = [1,2,5,10,20,50,100]

matchThreshold = -0.5

closeThreshold = 0.1

#http://my.safaribooksonline.com/book/programming/python/0596007973/python-shortcuts/pythoncook2-chp-4-sect-16
correctnessDict = {
	"species":	checkCorrectSpeciesOlsonFormat,
	"genus":	checkCorrectGenusOlsonFormat,
}

scoringMethodDict = {
	"raiphy":	scoreRAIphy,
}