#! /usr/bin/env python

import sys, os, getopt, time

def main(argv):
	myfile = ''
	out = ''
	try:
		opts, args = getopt.getopt(argv,"f:o:p:")
	except getopt.GetoptError:
		sys.exit(2)
		
	for opt, arg in opts:
		if opt == "-f": 
			myfile = arg
		elif opt == "-o":
			out = arg	
		elif opt == "-p":
			path = arg
		#elif opt == "-a":
		#	pref = arg
	
	outlog = out + "_log"
	fOut = open(outlog,'w')
	fOut.close()
	
	#jList = [0.5,0.7,0.9]
	#jList = [0.7,0.8,0.9]
	jList = [0.9]
	#nList = [0.01,0.03,0.1]
	#nList = [0.02,0.05,0.1]
	#nList = [0.08,0.1,0.12]
	nList = [0.1]
	#tDict = {'allClose': [2,4,6,8,10,12,14,16,18,20,22],
	#		 'lowClose': [2,4,6,8,10,14,18,22],
	#		 '2by4': [2,6,10,14,18,22],
	#		 '4by4': [4,8,12,16,20]}
	#tDict = {'4allClose': [4,6,8,10,12,14,16,18,20,22],
	#		 '4lowClose': [4,6,8,10,14,18,22],
	#		 '4by4': [4,8,12,16,20]}
	#tDict = {'4allCloseChop': [4,6,8,10,12,14,16,18],
	#		'4lowCloseChop': [4,6,8,10,14,18]}
	#tDict = {'4by4Chop': [4,8,12,16]}
	tDict = {'4allCloseChop': [4,6,8,10,12,14,16,18]}
	#sDict = {'tacoa': [0.3,0.35,0.4],
	#		'tetra': [0.2,0.4,0.6],
	#		'raiphy': [-17.0,-16.5,-16.0]}
	#sDict = {'tetra': [0.6,0.7],
	#		'raiphy': [-17.0,-16.0]}
	sDict = {'tetra': [0.6,0.7]}
	#prefList = ['min','mean','median','max']
	#prefList = ['max','90','80','70','len']
	prefList = ['60','median','40','len2']
	#for score in sorted(sDict.keys()):
	
	try:
		for pref in prefList:
			for t in tDict:			
				for score in sDict.keys():
					for s in sDict[score]:
						for j in jList:
							for n in nList:
		#for pref in ['median']:
		#	for t in ['4allCloseChop']:
		#		for s in [-17.0]:
		#			for j in [0.7]:
		#				for n in [0.1]:
								splitList = [s]*(len(tDict[t])-1)
								tStr = "["+",".join([str(ti) for ti in tDict[t]])+"]"
								sStr = "["+",".join([str(si) for si in splitList])+"]"
								print "{!s}, {!s}, {!s}, {!s}, {!s}, {!s}". format(score, n, j, tStr, sStr, pref)
								myOut = "{!s}_{!s}_N_{!s}_J_{!s}_C_{!s}_L_{!s}_A_{!s}".format(out,score,n,j,t,s,pref)
								os.system('echo "::{!s},{!s},{!s},{!s},{!s},{!s}" >> {!s}'.\
									format(score, n, j, tStr, sStr, pref, outlog))
								start = time.time()
								os.system("python horatio.py -i {!s} -o {!s} -s {!s} -n {!s} -j {!s} -c {!s} -l {!s} -p {!s} -a {!s} >> {!s}".\
									format(myfile, myOut, score, n, j, tStr, sStr, path, pref, outlog))
								end = time.time()
								os.system('echo "{:03.2f}" >> {!s}'.format(end-start, outlog))
	except KeyboardInterrupt:
		sys.exit()
	#fOut.close()			

if __name__ == "__main__":
	main(sys.argv[1:])