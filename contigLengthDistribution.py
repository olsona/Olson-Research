from matplotlib import pyplot
import math

def dist(infile):
	f = open(infile,'r')
	lines = f.readlines()
	num = len(lines)
	f.close()
	lens = [0.0 for i in range(num/2)]
	max = 0
	for i in range(num/2):
		lens[i] = len(lines[(i*2)+1])
		if lens[i] > max:
			max = lens[i]
			print max
	mmax = int(math.ceil(max/1000))
	bins = [i*1000 for i in range(mmax+1)]
	pyplot.hist(lens,bins,log=True)
	pyplot.show()