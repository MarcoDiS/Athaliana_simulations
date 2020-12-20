import sys
import matplotlib
matplotlib.use("Agg")

from matplotlib import pyplot as plt
import numpy as np
from pytadbit.utils.extraviews import nicer
from pytadbit.parsers.hic_parser import read_matrix

resolution = 1
matrix = read_matrix("raw_matrix.txt")

#print len(matrix),"x",len(matrix[0])
#print matrix
        
print "Saving the raw matrix"
from pytadbit.mapping.analyze import hic_map
print "I am saving the hic map both in a text file and in a figure" 
sys.stdout.flush()

filename="raw_matrix.txt"
fp = open(filename, "r")

size = int(sys.argv[1])

#matrix = np.zeros((size,size))
#i=0
#for line in fp.readlines():
#    line  = line.strip()
#    splitt = line.split()
#    for j in xrange(len(splitt)):
#        matrix[i][j] = float(splitt[j])
#    i+=1

#plt.imshow(matrix, interpolation="None", origin="lower") #, vmin=0, vmax=1.0)
#plt.title("Contact map at %s" % (nicer(resolution)))                                                                      
#plt.colorbar(label='Contacts')                                                                     
#plt.savefig('raw_matrix.pdf')
#plt.close("all")

hic_map(matrix, 
        resolution = resolution,
        #clim       = (-10,10),
        savefig    = "raw_matrix.pdf",
        normalized = False)
    
