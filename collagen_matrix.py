#!/usr/bin/env python

'''
Created on May 3rd, 2017
author: Jennifer Teubl

This script takes as input a file containing the results from an in silico trypsin digest.
Species scientific names have their own line starting with '>' and followed by each peptide that can be mapped to it.
A pairwise comparison between all species is done by checking for the number of common peptides two species have,
the result is used as a metric for a heatmap (computed in R)

input 1: peptide sequence text file sorted by species
input 2: a text file with scientific names and common names. Columns are tab separated.
input 3: path and file name for the output matrix

output1: matrix containing pairwise calculations
'''

import sys
import pandas as pd
import itertools
import numpy as np
import os
import re


seqFile = sys.argv[1] #input 1
names = sys.argv[2] #input 2
output = sys.argv[3] #input 3

n = open(names, 'rb')
nameDict = {}
for line in n:
    sci, com = line.strip().split('\t')
    nameDict[sci] = com

def openFile(f):
    speciesPeptides = {}
    peptides = False
    firstline = False
    for line in f:
        if line.startswith('>'):
            name = line[1:].strip()
            if peptides:
                speciesPeptides[firstline]=peptides
                peptides = []
            else:
                peptides = []
        else:
            peptides.append(line.strip())
            firstline = name
    speciesPeptides[firstline] = peptides
    for i in speciesPeptides:
        vals = speciesPeptides[i]
        vals = sorted(vals)
        vals = list(set(vals))
        speciesPeptides[i] = vals
    
    return speciesPeptides

f = open(seqFile, 'rb')
speciesPeptides = openFile(f)


speciesNames = speciesPeptides.keys()
rowNames = []
colnames = []
for species in speciesNames:
    numPepPos = len(speciesPeptides[species])
    com = nameDict[species]
    colnames.append(com)
    if numPepPos <100:
        rowNames.append("%i\t\t\t\t  %s (%s)"%(numPepPos, com, species))
    else:
        rowNames.append("%i\t\t\t\t%s (%s)"%(numPepPos, com, species))

similarityRatio = []
for keyTest in speciesNames:
    test = speciesPeptides[keyTest]
    ratios = []
    for keyComp in speciesNames:
        comp = speciesPeptides[keyComp]
        total = len(comp) + len(test)
        dup = 2 * len(list(set(comp) & set(test)))
        ratio = float(dup)/float(total)
        ratios.append("%.2f"%round(ratio,2))
    similarityRatio.append(ratios)
    

dfRatio = pd.DataFrame(similarityRatio, columns = colnames, index = rowNames)
with open(output, 'wb') as out:
    dfRatio.to_csv(out, sep = '\t')

