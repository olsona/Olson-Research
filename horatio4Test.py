#! /usr/bin/env python

import sys, os, getopt, time

def main(argv):
	myfile = ''
	out = ''
	try:
		opts, args = getopt.getopt(argv,"f:o:p:s:")
	except getopt.GetoptError:
		sys.exit(2)
		
	for opt, arg in opts:
		if opt == "-f": 
			myfile = arg
		elif opt == "-o":
			out = arg	
		elif opt == "-p":
			path = arg
		elif opt == "-s":
			score = arg
	
	outlog = out + "_" + score + "_log"
	fOut = open(outlog,'w')
	
	jList = [0.5,0.7,0.9]
	#nList = [0.01,0.03,0.1]
	nList = [0.02,0.05,0.1]
	#tDict = {'allClose': [2,4,6,8,10,12,14,16,18,20,22],
	#		 'lowClose': [2,4,6,8,10,14,18,22],
	#		 '2by4': [2,6,10,14,18,22],
	#		 '4by4': [4,8,12,16,20]}
	#tDict = {'4allClose': [4,6,8,10,12,14,16,18,20,22],
	#		 '4lowClose': [4,6,8,10,14,18,22],
	#		 '4by4': [4,8,12,16,20]}
	tDict = {'4allCloseChop': [4,6,8,10,12,14,16,18],
			'4lowCloseChop': [4,6,8,10,14,18]}
			#'4by4Chop': [4,8,12,16]}
	sDict = {'tacoa': [0.3,0.35,0.4],
			'tetra': [0.2,0.4,0.6],
			'raiphy': [-17.0,-16.5,-16.0]}
	prefList = ['min','mean','median','max']
	#for score in sorted(sDict.keys()):
	sList = sDict[score]
	
	
	for t in tDict:			
		for s in sList:
			for j in jList:
				for n in nList:
					for pref in prefList:
		#for t in ['4by4Chop']:
		#	for s in [-17.0]:
		#		for j in [0.5]:
		#			for n in [0.02]:
						try:
							splitList = [s]*(len(tDict[t])-1)
							tStr = "["+",".join([str(ti) for ti in tDict[t]])+"]"
							sStr = "["+",".join([str(si) for si in splitList])+"]"
							print "{!s}, {!s}, {!s}, {!s}, {!s}, {!s}". format(score, n, j, tStr, sStr, pref)
							myOut = "{!s}_{!s}_N_{!s}_J_{!s}_C_{!s}_L_{!s}_A_{!s}".format(out,score,n,j,t,s,pref)
							fOut.write("{!s},{!s},{!s},{!s},{!s},{!s}".\
								format(score, n, j, tStr, sStr, pref))
							start = time.time()
							os.system("python horatio.py -i {!s} -o {!s} -s {!s} -n {!s} -j {!s} -c {!s} -l {!s} -p {!s} -a {!s}".\
								format(myfile, myOut, score, n, j, tStr, sStr, path, pref))
							end = time.time()
							fOut.write("{:03.2f}\n".format(end-start))
							minutes = int((end-start)/60.0)
							secs = int((end-start)-60.0*minutes)
							print "{!s}:{!s}\n".format(minutes,str(secs).zfill(2))
						except KeyboardInterrupt:
							sys.exit()
	fOut.close()			

if __name__ == "__main__":
	main(sys.argv[1:])