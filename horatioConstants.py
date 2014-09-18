import horatioFunctions as hfun
import horatioUtils as hutil
import numpy

myMethods = ['raiphy']

usageString = 'bootstrap.py -i <input file> -o <output file> -p <score computation path> '\
			'[-c <cut schedule> -s <scoring function>]'
				
bigUsageString = '-i, --ifile\t\tInput metagenomics FASTA file\n'+\
				'-o, --ofile\t\tOutput file\n'+\
				'-c, --cut\t\tPartition schedule (optional)\n'+\
				'-s, --score\t\tScoring function for use (optional)\n'+\
				'-p, --path\t\tPath to extra code for computation (only required if scoring function is RAIphy)\n'+\
				'\t\t\t	 format: [1,2,5]\n'+\
				'\t\t\t	 [1,2,5] implies the partition schedule > 5kbp, 2kbp, 1kbp\n'
				
defaultSchedule = [4,6,8,10,12,14,16,18]

#http://my.safaribooksonline.com/book/programming/python/0596007973/python-shortcuts/pythoncook2-chp-4-sect-16
# ***
#correctnessDictSeed = {
#	 "species": checkCorrectSpeciesOlsonFormat,
#	 "genus":	checkCorrectGenusOlsonFormat,
#}

#correctnessDictMax = {
#	 "species": checkCorrectMatchClusterMaxSpecies,
#	 "genus":	checkCorrectMatchClusterMaxGenus,
#}
# ***

scoringMethod = {
	"raiphy":	hfun.scoreRAIphy,
	"tetra":	hfun.scoreTETRA,
	"tacoa":	hfun.scoreTACOA,
}

scoringMethodFinal = {
	"raiphy":	hfun.scoreRAIphyFinal,
	"tetra":	hfun.scoreTETRAFinal,
	"tacoa":	hfun.scoreTACOAFinal,
}

apPreferences = {
	"min":	numpy.min,
	"median":	numpy.median,
	"mean":	numpy.mean,
	"max":	numpy.max,
	"90":	hutil.percentile90,
	"80":	hutil.percentile80,
	"70":	hutil.percentile70
}