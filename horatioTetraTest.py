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
    nList = [0.01,0.02,0.05,0.1,0.2]
    tList = [[2,4,6,8,10,12,14,16,18,20,22],\
                [2,4,6,8,10,14,18,22],\
                [2,6,10,14,18,22],\
                [4,8,12,16,20]]
    sList = [0.4,0.45,0.5,0.55,0.6]
    for j in jList:
        print j  
        for n in nList:
            for t in tList:
                for s in sList:
                    splitList = [s]*len(t)
                    #print "python horatio.py -i {!s} -o {!s} -c {!s}"\
                    #    +"-s tetra -n {!s} -j {!s} -l {!s}".format(myfile,\
                    #    out, t, n, j, splitList)
                    #os.system("python horatio.py -i {!s} -o {!s} -c {!s}"\
                    #    +"-s tetra -n {!s} -j {!s} -l {!s}".format(myfile,\
                    #    out, t, n, j, splitList))

if __name__ == "__main__":
    main(sys.argv[1:])