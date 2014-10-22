#! /usr/bin/env python

import sys, os, getopt, time
import os.path

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
			scr = arg
	
	outlog = out + "_log"
	fOut = open(outlog,'a')
	fOut.close()
	
	jList = [0.5,0.7,0.9]
	#jList = [0.7,0.8,0.9]
	
	#nList = [0.01,0.03,0.1]
	#nList = [0.01,0.05,0.1]
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
	#		'4lowCloseChop': [4,6,8,10,14,18],
	#		'4by4Chop': [4,8,12,16]}
	#tDict = {'4lowCloseChop': [4,6,8,10,14,18]}
	tDict = {'4allCloseChop': [4,6,8,10,12,14,16,18]}
	
	#sDict = {'tacoa': [0.3,0.35,0.4],
	#		'tetra': [0.2,0.4,0.6],
	#		'raiphy': [-17.0,-16.5,-16.0]}
	scoreDict = {'RA':'raiphy','TA':'tacoa','TE':'tetra'}
	score = scoreDict[scr]
	#sDict = {'tetra': [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7],
	#		'raiphy': [-18.0,-17.5,-17.0,-16.5,-16.0,-15.5,-15.0],
	#		'tacoa': [0.2,0.24,0.28,0.32,0.36,0.4,0.44]}
	sDict = {'tetra': [0.5,0.55,0.6]}
	#sDict = {'tetra': [0.6],
	#		'raiphy': [-17.0],
	#		'tacoa': [0.4]}
			
	#prefList = ['min','mean','median','max']
	#prefList = ['max','90','80','70','len']
	#prefList = ['60','median','40','len2']
	#prefList = ['median']
	prefList = ['min','mean','median','max','95','90','80','70','60','40','len','len2','len3','len4']
	#prefList = ['40','60','max']
	#prefList = ['max']
	#prefList = ['len4']
	#prefList = ['95','len3']
	#prefList = ['N']
	#for score in sorted(sDict.keys()):
	score='tetra'
	
	try:
		for pref in prefList:
			for t in tDict:			
				for s in sDict[score]:
					for j in jList:
						for n in nList:
							splitList = [s]*(len(tDict[t])-1)
							tStr = "["+",".join([str(ti) for ti in tDict[t]])+"]"
							sStr = "["+",".join([str(si) for si in splitList])+"]"
							myOut = "{!s}_{!s}_N_{!s}_J_{:0.2f}_C_{!s}_L_{:0.2f}_A_{!s}".format(out,score,n,j,t,s,pref)
							if(not(os.path.isfile(myOut+"_clusters_pickle"))):
								print "{!s}, {!s}, {!s}, {!s}, {!s}, {!s}". format(score, n, j, tStr, sStr, pref)
								os.system('echo "{!s},{!s},{!s},{!s},{!s},{!s}" >> {!s}'.format(score, n, j, t, s, pref, outlog))
								start = time.time()
								#os.system("python horatio.py -i {!s} -o {!s} -s {!s} -n {!s} -j {!s} -c {!s} -l {!s} -d 1 -p {!s} -k 2 -a {!s} >> {!s}".format(myfile, myOut, score, n, j, tStr, sStr, path, pref, outlog))
								os.system("python horatio.py -i {!s} -o {!s} -s {!s} -n {!s} -j {!s} -c {!s} -l {!s} -d 0 -p {!s} -k 2  >> {!s}".format(myfile, myOut, score, n, j, tStr, sStr, path, outlog))
								end = time.time()
								os.system('echo "{:03.2f}::" >> {!s}'.format(end-start, outlog))
	except KeyboardInterrupt:
		sys.exit()
	#fOut.close()			

if __name__ == "__main__":
	main(sys.argv[1:])