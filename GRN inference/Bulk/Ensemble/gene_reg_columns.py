# -*- coding: utf-8 -*-
"""gene-reg-columns.ipynb

Script for input geneclust_jaccard.pl

"""

import sys
import re
import os
import pandas as pd

clust = {}
file_input ="ensemble_netw_MDD.txt"
file_output = "gene_reg_colMDD.txt"
f = open(file_input, 'r')
for line in f.readlines():
	fields = line.strip('\n').split('\t')
	gene = fields[1]
	#gene = gene.strip('\"')
	regulator = fields[0]
	if gene in clust.keys():
		clust[gene].append(regulator)
	else:
		listfeature = []
		listfeature.append(regulator)
		clust[gene] = listfeature

f_out = open(file_output, 'w')
genelist = pd.DataFrame.from_dict(clust, orient ='index')
genelist = genelist.iloc[1: , :]
print(genelist.head())
genelist.to_csv(f_out, header = False, sep = '\t')