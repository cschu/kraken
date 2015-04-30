#!/usr/bin/env python

import os
import sys
import csv

from collections import Counter

from kraken_visualize import readTaxonomyNames, readTaxonomyNodes, getDescendents


taxons = readTaxonomyNames(os.path.join(sys.argv[2], 'names.dmp'))

taxID = int(sys.argv[3])
validTaxons = getDescendents(taxID, readTaxonomyNodes(os.path.join(sys.argv[2], 'nodes.dmp'))[1])

c = Counter([int(row[2]) 
	         for row in csv.reader(open(sys.argv[1]), delimiter='\t', quotechar='"')])

N = float(sum(c.values()))
ct = 0
for k in sorted(c, key=lambda x:c[x], reverse=True):
	if k in validTaxons:
		print k, taxons.get(k), c[k], '%.10f' % (c[k]/N)
		ct += c[k]
	# print k, taxons.get(k, 'N/A'), c[k], 'VALID' if k in validTaxons else ''
print ct, ct/N
