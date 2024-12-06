
"""
# @author: Hanne Puype
Created on Wed 26 July 2023
"""

#### 
# Script to use unibind as ground truth for networks
# Unibind data comes as TFBS, so we need to convert it to genes 
# Downloaded Unibind data from: http://unibind.uio.no/downloads Homo sapiens
# We can use pybiomart for this
# For the genes, there are two option: chose the closest gene or take a certain distance
# For the distance, we can use pybedtools

# After this, we can compare the networks with the unibind data with precision and recall 
#####

import os
import sklearn as sk
import numpy as np
import pandas as pd
import seaborn as sns
import pybedtools
import pybiomart as pbm

# path to Unibind file
Unibind_path = "/home/hanne/Downloads/"
Unibind_file = "hg38_TFBSs_TF.bed"


def prepare_TSS_bed_file(biotype_filter = 'protein_coding'):
    """
    Download and filter transcription start sites from emsembl biomart, 
    and return a bed file with TSS positions, as well as a dictionary of all genes in the file
    """
    # Get data from ensembl
    dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
    # Define attributes to receive
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    # Filter out unwanted entires
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
    filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    # replace non-uniform naming
    def add_chr_prefix_to_chrom(df):
        # Update 'chrom' column with 'chr' prefix if the column starts with a digit
        df.loc[df['Chromosome/scaffold name'].str.contains('^\d', regex=True), 'Chromosome/scaffold name'] = 'chr' + df['Chromosome/scaffold name']
        # update 'chrom' column with 'chrUn' prefix if column starts with 'KI'
        df.loc[df['Chromosome/scaffold name'].str.contains('^KI', regex=True), 'Chromosome/scaffold name'] = 'chrUn_' + df['Chromosome/scaffold name']
        # change 'chrUn_KI270734.1' to 'chrUn_KI270734v1'
        df.loc[df['Chromosome/scaffold name'].str.contains('chrUn_KI270734.1', regex=True), 'Chromosome/scaffold name'] = 'chrUn_KI270734v1'
        # change 'X' to 'chrX' and 'Y' to 'chrY'
        df.loc[df['Chromosome/scaffold name'].str.contains('Y', regex=True), 'Chromosome/scaffold name'] = 'chrY'
        df.loc[df['Chromosome/scaffold name'].str.contains('X', regex=True), 'Chromosome/scaffold name'] = 'chrX'
    add_chr_prefix_to_chrom(annot)   
    # Set column names
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    # filter for biotype protein coding 
    if biotype_filter != None:
        annot = annot[annot.Transcript_type == biotype_filter]
    annot.insert(2, "End", annot['Start']) # add End column, which is the same as Start column
    annot.dropna(inplace=True)
    
    # get all non NaN gene names in dict
    genes = {}
    for i in annot['Gene']:
        if i not in genes and i != 'NaN':
            genes[i] = ''
    
    # Make bedfile, retunr it and also save result as bed file
    x = pybedtools.BedTool.from_dataframe(annot)

    return x, genes

# prepare TSS bed file
TSS_bed, genes = prepare_TSS_bed_file()
print("Number of genes:", len(genes))
print("Number of TSS:", len(TSS_bed))
print(TSS_bed.head()) # chr, start, end, strand, gene, transcript type
TSS_bed = TSS_bed.sort()
print(TSS_bed.head())

TSS_bed_df = TSS_bed.to_dataframe()
print(TSS_bed_df.head())
TSS_bed_df['chrom'].unique()
TSS_bed = pybedtools.BedTool.from_dataframe(TSS_bed_df)
TSS_bed_df.to_csv('/home/hanne/MyFiles/Master/TSS.bed', index=False, sep='\t', header=False)

# used bedtools closest on command line
# bedtools closest -d -a hg38_TFBSs_TF.bed -b TSS.bed > closest_UniBind_TSS.bed

# filter bed file for duplicates and only keep columns needed (TF, gene, distance) on command line
# awk -F'\t' '!seen[$4,$14]++ {print $4, $14, $16}' closest_UniBind_TSS.bed > closest_UniBind_TSS_ndup.bed

