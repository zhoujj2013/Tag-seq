import os, sys
import re
import numpy as np
import scipy.signal
#import six

def usage():
    print('\nCombine cufflink quantification result to matrix\n')
    print('Author: zhoujj2013@gmail.com 8/29/2016\n')
    print('Usage: python '+sys.argv[0]+' 1.gene_tracking:aa 2.gene_tracking:bb ...')
    print('')
    sys.exit(2)

# check args
if len(sys.argv) < 3:
    usage()


dep_plus_f = sys.argv[1]
dep_minus_f = sys.argv[2]

min_dep = int(sys.argv[3])
min_dis = int(sys.argv[4])

#prefix = sys.argv[5]

dep_plus = open(dep_plus_f, 'rb')
dep_minus = open(dep_minus_f, 'rb')

# plus strand
plus = {}
plus_dep = {}
while True:
    l = dep_plus.readline()
    if len(l) == 0:
        break
    lc = l.strip("\n").split("\t")
    nc = lc[3].split("_")
    name_index = nc.pop()
    id_str = "_".join(nc)
    if id_str in plus:
        plus[id_str].append(lc)
    else:
        plus[id_str] = []
        plus[id_str].append(lc)
 
    if id_str in plus_dep:
        plus_dep[id_str].append(lc[-1])
    else:
        plus_dep[id_str] = []
        plus_dep[id_str].append(lc[-1])

for k, v in plus_dep.items():
    vector = np.array(v)
    indexes, _ = scipy.signal.find_peaks(vector, height=min_dep, distance=min_dis)
    if len(indexes) == 0:
        continue
    else:
        for i in indexes:
            arr = plus[k]
            print "\t".join(arr[i])
# minus strand
minus = {}
minus_dep = {}
while True:
    l = dep_minus.readline()
    if len(l) == 0:
        break
    lc = l.strip("\n").split("\t")
    nc = lc[3].split("_")
    name_index = nc.pop()
    id_str = "_".join(nc)
    if id_str in minus:
        minus[id_str].append(lc)
    else:
        minus[id_str] = []
        minus[id_str].append(lc)
 
    if id_str in minus_dep:
        minus_dep[id_str].append(lc[-1])
    else:
        minus_dep[id_str] = []
        minus_dep[id_str].append(lc[-1])

for k, v in minus_dep.items():
    vector = np.array(v)
    indexes, _ = scipy.signal.find_peaks(vector, height=min_dep, distance=min_dis)
    if len(indexes) == 0:
        continue
    else:
        for i in indexes:
            arr = minus[k]
            print >>sys.stderr, "\t".join(arr[i])

#print plus_dep
#print >>sys.stderr, plus

