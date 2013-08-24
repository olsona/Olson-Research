def readSequence(fi):
    f = open(fi,'r')
    concat = ''
    buf = f.readline().rstrip()
    while buf:
        seq_name, seq = buf[2:], ''
        buf = f.readline()
        while not buf.startswith('>') and buf:
            seq = seq + buf.rstrip()
            buf = f.readline()
        if seq_name.find("complete genome") == -1:
            concat = concat + seq
        elif seq_name.find("complete genome") > -1:
            concat = seq
            break
    return seq_name, concat
    f.close()

def chopRandomTax(seq, min_size, num_chop, path, name, taxonomy):
    import random
    llim = len(seq) - min_size - 300
    if llim <= 0:
        return
    if min_size >= 1000:
        sizestr = str(min_size/1000) + "k"
    else:
        sizestr = str(min_size)
    f = open(path + "/" + name + "_random_chopped_" + sizestr + "b.fna",'w')
    for i in range(num_chop):
        st = random.randrange(llim)
        leng = random.randrange(300) + min_size
        subseq = seq[st:st+leng]
        f.write(">:" + taxonomy + ": " + str(st) + " - " + str(st + leng) + "\n")
        ct = 0
        while 1:
            f.write(subseq[ct:ct+60] + "\n")
            ct += 60
            if (ct > len(subseq)):
                break
        f.write("\n")
    f.close()
    
def chopRandomInOrder(seq, size, name, out):
    import random
    interval = size/10
    ct = 0
    file = open(out,'w')
    while ct < len(seq)+interval:
        mylen = size + random.randrange(-interval, interval)
        myseq = seq[ct:ct+mylen]
        file.write(">{!s}; {!s}-{!s}\n".format(name,ct,ct+mylen))
        file.write("{!s}\n".format(myseq))
        ct += mylen
    file.close()

def chopAllFolder(size, path):
    import os
    fls = os.listdir(path)
    for f in fls:
        if f[-3:] == 'fna':
            nm = f[:-4]
            tax, se = readSequence(path+"/"+f)
            chopRandomInOrder(se, size, nm, "{!s}/Chopped/{!s}_{!s}.fna".format(path, nm, size))
 

def writeProperSequences(path1, path2):
    import os
    fls = os.listdir(path1)
    for f in fls:
        if f[-3:] == 'fna':
            f2 = open(path2 + "/" + f, 'w')
            nm, se = readSequence(path1 + "/" + f)
            f2.write(">: " + nm + "\n")
            ct = 0
            while 1:
                f2.write(se[ct:ct+60] + "\n")
                ct += 60
                if (ct > len(se)):
                    break
            f1.close()
            f2.close()