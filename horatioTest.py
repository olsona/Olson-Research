#! /usr/bin/env python

import sys, os, getopt

def main(argv):
    myfile = ''
    out = ''
    try:
	opts, args = getopt.getopt(argv,"f:o:")
    except getopt.GetoptError:
	sys.exit(2)
        
    for opt, arg in opts:
        if opt == "-f":	
            myfile = arg
        elif opt == "-o":
            out = arg	
    
    outlog = out + "_log"
    fOut = open(outlog,'w')
    
    jList = [0.5,0.7,0.9]
    nList = [0.01,0.03,0.1]
    tDict = {'allClose': [2,4,6,8,10,12,14,16,18,20,22],
            'lowClose': [2,4,6,8,10,14,18,22],
            '2by4': [2,6,10,14,18,22],
            '4by4': [4,8,12,16,20]}
    sDict = {'tacoa': [0.3,0.4,0.5],
            'tetra': [0.2,0.4,0.6,0.8],
            'raiphy': [-17.0,-16.0,-15.0,-14.0]}
    for score in sDict:
        sList = sDict[score]
        for t in tDict:
            for s in sList:
                for j in jList:
                    for n in nList:
                        splitList = [s]*(len(tDict[t])-1)
                        tStr = "["+",".join([str(ti) for ti in tDict[t]])+"]"
                        sStr = "["+",".join([str(si) for si in splitList])+"]"
                        myOut = "{!s}_{!s}_N_{!s}_J_{!s}_C_{!s}_L_{!s}".format(out,score,n,j,t,s)
                        print "-s {!s} -n {!s} -j {!s} -c {!s} -l {!s}".\
                            format(score, n, j, tStr, sStr)
                        os.system("python horatio.py -i {!s} -o {!s} -s {!s} -n {!s} -j {!s} -c {!s} -l {!s}".\
                            format(myfile, myOut, score, n, j, tStr, sStr))

    fOut.close()            

if __name__ == "__main__":
    main(sys.argv[1:])