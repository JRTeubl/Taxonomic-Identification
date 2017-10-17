#!/usr/bin/env python

'''
Created on November 4th, 2016
author: Jennifer Teubl

This program quantifies and visualizes database metrics. Takes as input a text file with protein digest and a protein name : species
file
'''

import sys
import matplotlib.pyplot as plt

import operator
import numpy as np
import itertools

seqFile = sys.argv[1]
output = sys.argv[2]
colType = sys.argv[3]

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
        speciesPeptides[i] = vals
    
    return speciesPeptides

f = open(seqFile, 'rb')
speciesPeptides = openFile(f)

def NumPepPerSpecies(speciesPeptides):
    #total number of non-redundant peptides per species
    speciesNumPep = {}
    for key in speciesPeptides:
        numPep = len(speciesPeptides[key])
        speciesNumPep[key] = numPep
    avg = np.mean(speciesNumPep.values())
    std = np.std(speciesNumPep.values())
    numPep = sorted(speciesNumPep.items(), key=operator.itemgetter(1))
    return numPep, avg, std

totalPepPerSpecies, totalAvg, totalStd = NumPepPerSpecies(speciesPeptides)

def ProteotypicPepPerSpecies(speciesPeties):
    # number of species proteotypic peptides (peptides not found in any other species)
    proteotypic = {}
    for key in speciesPeptides:
        vals = speciesPeptides[key]
        compare = []
        for keyNew in speciesPeptides:
            if key != keyNew:
                compare.append(speciesPeptides[keyNew])
        compare = list(set(itertools.chain.from_iterable(compare)))
        keep = []
        for x in vals:
            if x not in compare:
                keep.append(x)
        proteotypic[key] = keep
    return proteotypic

proteotypic = ProteotypicPepPerSpecies(speciesPeptides)
proteoPerSpecies, proteoAvg, proteoStd = NumPepPerSpecies(proteotypic)
totalPepPerSpecies = dict(totalPepPerSpecies)

#sort bar plots by the number of species proteotypic peptides
fig, ax1 = plt.subplots(1, figsize=(10,8))
proteo = [i[1] for i in proteoPerSpecies]
names = [i[0] for i in proteoPerSpecies]
total = [totalPepPerSpecies[i[0]] for i in proteoPerSpecies]

bar_width = 1
bar_l = [i for i in range(len(total))]
tick_pos = [i+.5 for i in bar_l]
ax1.bar(bar_l, total, bar_width, color = 'DarkOrange', alpha = .7, label = 'Total Peptide Count')
ax1.bar(bar_l, proteo, bar_width, color = 'DarkOliveGreen', alpha = 1, label = 'Proteotypic Peptide Count')
#, bottom = total

plt.xticks(tick_pos, names, rotation = 90)
plt.legend(loc='upper center', fontsize = 12)
plt.title('Peptides Per Species:%s'%colType)
plt.xlim([min(tick_pos)-bar_width, max(tick_pos)+bar_width])
ax1.set_ylabel("Number of Peptides")
ax1.set_xlabel("Species")
fig.tight_layout()
# plt.show()
plt.savefig(output+colType+'_peptideMatrics.pdf')


