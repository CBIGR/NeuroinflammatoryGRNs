# -*- coding: utf-8 -*-
"""
Convert_ModuleViewer
Right format for module genes and regulators to load into ModuleViewer

"""

import sys
import re
import os

clust = {}
file_input = "reg_clust_ens_AD.txt"
file_output = "kmed_reg_MV_AD.txt"
f = open(file_input, 'r')
for line in f.readlines():
	fields = line.strip('\n').split('\t')
	gene = fields[0]
	gene = gene.strip('\"')
	cluster = fields[1]
	if cluster in clust.keys():
		clust[cluster].append(gene)
	else:
		listfeature = []
		listfeature.append(gene)
		clust[cluster] = listfeature
f_out = open(file_output, 'w')

for c in clust.keys():
	ids = clust[c]
	idlist = '|'.join(ids)
	line = ''.join([c,"\t",idlist])
	f_out.write(line+"\n")