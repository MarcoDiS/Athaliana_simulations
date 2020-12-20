import numpy
import matplotlib
import os
import scipy
matplotlib.use("Agg")

import pylab as pl
from numpy import zeros, ones, arange, asarray, concatenate, dot

import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
import scipy as sp
from scipy.optimize import linprog
from math import sqrt

from scipy.spatial import ConvexHull
import sys


def pnt_in_cvex_hull(hull, point, tolerance=1e-3):
    return all(
        (dot(eq[:-1], point) + eq[-1] <= tolerance)
        for eq in hull.equations)

threshold = float(sys.argv[1]) # overlap threshold
npoints = int(sys.argv[2])

radius=0 

cluster = []
for i in xrange(npoints):
    cluster.append(i)

infile=open("_overlaps","r")
distances = zeros((npoints,npoints),dtype=float)
for line in infile.readlines():
    line = line.strip().split()
    distances[int(line[0])][int(line[1])] = float(line[2])
    distances[int(line[1])][int(line[0])] = float(line[2])


for point1 in xrange(npoints):
    for point2 in xrange(point1+1,npoints):
        if cluster[point2] == cluster[point1]:
            #print "Same cluster of %s, skipping %s" % (point1, point2)
            continue
        if distances[point1][point2] > (threshold+2.0*radius):
            
            # Take the lower cluster index
            c = cluster[point1]
            if cluster[point2] < cluster[point1]:
                c = cluster[point2]

            # Join the clusters of point1 and point2 assigning the lower cluster index to all the elements
            for i in xrange(npoints):
                if cluster[i] == cluster[point2] or cluster[i] == cluster[point1]:
                    cluster[i] = c

            #print cluster[point1], cluster[point2]

#print cluster
cluster_set = set(cluster)
#print cluster_set
print len(cluster_set)

#max_npts=0
# Volume per cluster
#for c in cluster_set:
#    npts=0
#    for i in xrange(len(cluster)):
#        if cluster[i] == c:
#            npts+=1
#    if npts == 0:
#        continue
#    print "Npoints",c,npts
