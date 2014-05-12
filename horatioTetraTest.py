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
    
    jList = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    #jList = [0.5]
    nList = [0.01,0.02,0.05,0.1,0.2]
    #nList = [0.02]
    #tList = [[2,4,6,8,10,12,14,16,18,20,22],\
                [2,4,6,8,10,14,18,22],\
    #            [2,6,10,14,18,22],\
    #            [4,8,12,16,20]]
    tList = [[4,8,12,16,20]]
    #sList = [0.4,0.45,0.5,0.55,0.6]
    sList = [0.5]
    for j in jList:
        for n in nList:
            for t in tList:
                for s in sList:
                    splitList = [s]*(len(t)-1)
                    myOut = "{!s}_N_{!s}_J_{!s}_C_{!s}_L_{!s}".format(out,n,j,t,splitList)
                    tStr = "["+",".join([str(ti) for ti in t])+"]"
                    sStr = "["+",".join([str(si) for si in splitList])+"]"
                    print "python horatio.py -i {!s} -o {!s} -s tetra -n {!s} -j {!s} -c {!s} -l {!s}".format(myfile,\
                        myOut, n, j, tStr, sStr)
                    os.system("time python horatio.py -i {!s} -o {!s} -s tetra -n {!s} -j {!s} -c {!s} -l {!s}".format(myfile,\
                        myOut, n, j, tStr, sStr))

if __name__ == "__main__":
    main(sys.argv[1:])