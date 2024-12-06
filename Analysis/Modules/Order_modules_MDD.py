#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 09:51:07 2022

@author: boris
"""

########################################################################################################################################
# After visualizing your LemonTree modules with ModuleViewer you end up with a lot of pdf's that all have to be manually inspected
# We want to automate this and prioritze the modules based on a good disease/control separation
# I will look at genes involved in each module and try to find those modules with the biggest difference in average expression
########################################################################################################################################

import os
import pandas as pd
from scipy.stats import mannwhitneyu
import shutil
import sys
#import numpy as np


work_dir = "/home/hanne/MyFiles/Master"
expression_location = "/home/hanne/MyFiles/Master/Processed_counts_MDD_GS.txt"
mvf = "/home/hanne/MyFiles/Master/Ens_Anno_Diagn_MDD" # module viewer file for annotation
plot_folder = "/home/hanne/MyFiles/Master/ModulesMDD"

# Read clusters file, create dictionary with cluster names and genes in the cluster
# The file is the format that was made for ModuleViewer
with open("/home/hanne/MyFiles/Master/kmed_modules_MV_MDD.txt", 'r') as handle:
    modules = {}
    for line in handle:
        line = line.rstrip()
        module = line.split()[0]
        genes_in_module = line.split()[1].split('|')
        modules[module] = genes_in_module
    print('Succesfully created modules dictionary')
    
# Read mvf file and create a dictionary that links sample names to case/control annotation
with open(mvf, 'r') as handle2:
    annotations = {}
    for line in handle2:
        if line.startswith('::'):
            pass
        else:
            samples = line[0:].rstrip().split('|')
    for element in samples:
        sample = element.split(':')
        annotations[sample[0]] = sample[1]
    print('Succesfully read mvf file information')
    
module_p_values = {}    
expression_data = pd.read_csv(expression_location, sep='\t')
print(expression_data.head())
print(annotations)

for module in modules: # Loop through each module
    genes = modules[module]
    data_to_use = expression_data[expression_data["GeneSymbol"].isin(genes)]
    data_to_use = data_to_use.iloc[:, 1:len(data_to_use.columns)] # column with Genesymbol gone
    means = data_to_use.mean(axis=0) # Calculate the average expression of all relevant genes per sample: axis = 0 means along the column
    # #print(data_to_use)
    #print(means)
    means = means.to_frame()
    means.columns = ['expression']
    means.index.name = 'sample'
    means.reset_index(inplace=True)
    
    group1_samples = [key for key, value in annotations.items() if value != '#117733'] # in mvf: if sample has red annotation --> one group, so patients
    group2_samples = [key for key, value in annotations.items() if value == '#117733']
    
    group1_expression = means[means['sample'].isin(group1_samples)]
    group2_expression = means[means['sample'].isin(group2_samples)]
    
    group1 = group1_expression['expression'].tolist()
    group2 = group2_expression['expression'].tolist()
    
    stat, p = mannwhitneyu(group1, group2) # Non-parametric alternative for t test
    #print(module, p)
    module_p_values[module] = p
    #break
    
# Now we need to rearrange the original module files according to ascending p-values
modules_ordered = {k: v for k, v in sorted(module_p_values.items(), key=lambda item: item[1])}
count = 0

if not os.path.exists(plot_folder + 'ordered'):
    os.makedirs(plot_folder + 'ordered')
    

for module in modules_ordered:
    old_filename = plot_folder + '/module_' + module + '.png'
    #print(old_filename)
    new_filename = plot_folder + 'ordered/' + str(count) + '_module' + module + '.png'
    #print(new_filename)
    count += 1
    shutil.copy(old_filename, new_filename)
    
    
print(f'Renamed and prioritized {count} gene modules')

# extract p-values from dictionary and save to file for every module
print(module_p_values.keys())
df = pd.DataFrame.from_dict(module_p_values, orient='index') # from dict to dataframe
df.to_csv('/home/hanne/MyFiles/Master/p_values_MDD.csv', header=True)

# FDR correction with Benjamini-Hochberg method
from statsmodels.stats.multitest import multipletests
p_values = list(module_p_values.values())
reject, p_corrected, a1, a2 = multipletests(p_values, alpha=0.1, method='fdr_bh')
print(reject)
print(p_corrected)
# from list to dict with module names as keys
module_p_corrected = {}
for i in range(len(reject)):
    module_p_corrected[list(module_p_values.keys())[i]] = p_corrected[i]
print(module_p_corrected)

df = pd.DataFrame.from_dict(module_p_corrected, orient='index') # from dict to dataframe
df.to_csv('/home/hanne/MyFiles/Master/p_values_FDR_MDD.csv', header=True)

# significant modules
significant_modules = {}
for key, value in module_p_corrected.items():
    if value < 0.1:
        significant_modules[key] = value
print(significant_modules)

