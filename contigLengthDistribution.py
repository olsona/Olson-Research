from matplotlib import pyplot

def dist(infile):
    f = open(infile,'r')
    lines = f.readlines()
    num = len(lines)
    f.close()
    lens = [0.0 for i in range(num/2)]
    for i in range(num/2):
        lens[i] = len(lines[i*2])
    pyplot.hist(lens)
    pyplot.show()
