from scipy.stats import spearmanr as corr
import sys

replica    = int(sys.argv[1])

a = []
fp = open("_exp_%s" % (replica), "r")
for line in fp.readlines():
    line = line.strip().split()
    if len(line) > 1 and abs(float(line[0])-float(line[1]))>1:
        a.append(float(line[2]))
    else:
        a.append(float(line[0]))

b = []
fp = open("_model_%s" % (replica), "r")
for line in fp.readlines():
    line = line.strip().split()
    if len(line) > 1 and abs(float(line[0])-float(line[1]))>1:
        b.append(float(line[2]))
    else:
        b.append(float(line[0]))

print corr(a,b)[0] , len(a) #, len(b)
