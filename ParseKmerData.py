def getKmerData(k, path, out):
	import math, sys, os, numpy, itertools
	
	allDivs = []
	allReps = []
	allSDvs = []
	allPUnu = []
	allNUnu = []
	
	files = os.listdir(path)
	outfile = open(out,'w')
	if path[-1] != "/":
		path = path + "/"
	if k <= 10:
		alphabet = ["a","c","g","t"]
		for f in files:
			ks = "%02d" % k
			if len(f) > 12 and f[-11:] == "_" + ks + "_out.txt":
				total, div, rep, res = processFileFull(path + f, k)
				if total > 0.0:
					allDivs.append(div)
					allReps.append(rep)
					arr = numpy.array(res)
					stdv = arr.std(dtype=numpy.float)/total
					mean = arr.mean(dtype=numpy.float)
					allSDvs.append(stdv)
					myPUnu = []
					myNUnu = []
					for i in range(len(res)):
						if res[i] > numpy.float(2.0)*stdv + mean:
							myPUnu.append(unHash(i,alphabet,k))
						elif res[i] < numpy.float(-1.0)*stdv + mean:
							myNUnu.append(unHash(i,alphabet,k))
					allPUnu.append(myPUnu)
					allNUnu.append(myNUnu)		
		outfile.write("Divergences:\n%s\n\n" %(str(allDivs)))
		outfile.write("Representation:\n%s\n\n" %(str(allReps)))
#		for kay in range(2,6):
#			myunion = set()
#			for x in itertools.combinations(range(len(allDivs)), kay):
#				subs = set(allPUnu[x[0]])
#				for y in range(1,kay):
#					subs = subs & set(allPUnu[x[y]])
#				myunion = myunion | subs
#			outfile.write("Union of all intersections of %d sets has pct representation = %f\n\n" %(kay, 100.0*float(len(myunion))/float(pow(4,k))))
		outfile.write("Standard Deviations / Length:\n%s\n\n" %(str(allSDvs)))
		outfile.write("Frequent Kmers (freq > 2 sigma):\n%s\n\n" %(str(allPUnu)))
		outfile.write("Infrequent Kmers (freq < -sigma):\n%s\n\n" %(str(allNUnu)))
	else:
		for f in files:
			ks = "%02d" % k
			if f[-11:] == "_" + ks + "_out.txt":
				div, rep = processFileQuick(path + f)
				allDivs.append(div)
				allReps.append(rep)
		outfile.write("Divergences:\n%s\n\n" %(str(allDivs)))
		outfile.write("Representation:\n%s\n\n" %(str(allReps)))
	outfile.close()

def listLex(alphabet, n):
	import itertools as it
	res = []
	for p in it.product(alphabet, repeat=n):
		res.append(''.join(p))
	return res
	
def comprehension_flatten(iter_lst):
	return list(item for iter_ in iter_lst for item in iter_)

def getHash(string, alphabet, k):
	import math
	ans = 0
	try:
		for i in range(len(string)):
			c = string[i]
			n = alphabet.index(c)
			ans += n * (len(alphabet)**(k-i-1))
	except ValueError:
		ans = -1
	return int(ans)
	
def unHash(hsh, alphabet, k):
	import math
	ans = ''
	if hsh != -1:
		for i in range(k):
			res = hsh%len(alphabet)
			hsh = hsh/len(alphabet)
			ans = alphabet[res] + ans
	return ans

def processFileFull(file, k):
	import math, numpy
	f = open(file, 'r')
	total = numpy.float(0.0)
	div = numpy.float(0.0)
	rep = numpy.float(0.0)
	res = [numpy.float(0.0)]*int(pow(4,k))
	l = f.readline()
	while l != '':
		lres = l.rstrip().split('\t')
		if len(lres) == 2:
			if lres[0] == 'Total:':
				total = numpy.float(lres[1])
			elif lres[0] == 'Kullback-Leibler divergence:':
				div = numpy.float(lres[1])
			elif lres[0] == 'Fraction of kmers that are represented:':
				rep = numpy.float(lres[1])
			else:
				ind = getHash(lres[0],['a','c','g','t'],k)
				if ind > -1:
					res[ind] = numpy.float(lres[1])
		l = f.readline()
	f.close()
	return total, div, min(rep,1.0), res

def processFileQuick(file):
	import math
	f = open(file, 'r')
	div = 0.0
	rep = 0.0
	l = f.readline()
	while l[0:4] == 'Read':
		l = f.readline()
	l = f.readline()
	lres = l.rstrip().split('\t')
	div = float(lres[1])
	l = f.readline()
	lres = l.rstrip().split('\t')
	rep = float(lres[1])
	f.close()
	return div, min(rep,1.0)
	
def getKmerSpectra(path,setOfSpecies,k,out):
	import math, numpy
	if path[-1] != "/":
		path = path + "/"
	ks = "_%02d_" % k
	f = []
	for s in setOfSpecies:
		myf = open(path+s+ks+"out.txt", 'r')
		f.append(myf)
	o = open(out,'w')
	dict = {}
	for mf in f:
		total = 1.0
		l = mf.readline()
		while l != '':
			lres = l.rstrip().split('\t')
			if len(lres) == 2:
				if lres[0] == 'Total:':
					total = numpy.float(lres[1])
				else:
					ind = getHash(lres[0],['a','c','g','t'],k)
					if ind > -1:
						if lres[0] in dict:
							dict[lres[0]].append(math.pow(4,k)*numpy.float(lres[1])/total)
						else:
							dict[lres[0]] = [math.pow(4,k)*numpy.float(lres[1])/total]
			l = mf.readline()
		mf.close()
	skeys = sorted(dict.keys())
	exp = numpy.float(100.0/(math.pow(4,k)))
	for sk in skeys:
		arr = numpy.array(dict[sk]+[numpy.float(0.0)]*(len(setOfSpecies)-len(dict[sk])))
		mean = arr.mean(dtype=numpy.float)
		stdv = arr.std(dtype=numpy.float)
		o.write(sk+"\tMean: %f\tStdDv: %f\n"%(mean,stdv))
	o.close()
	