'''
This script takes as input a file linking species name with protein name and a digest file with peptides listed under protein names.
The output of this script is a file with peptides listed under species names (combined from all proteins)
'''

import os
import sys
import itertools

speciesFile = sys.argv[1]
peptideFile = sys.argv[2]
output = sys.argv[3]

protSpecies = {}
f = open(speciesFile, 'rb')
for line in f:
    p,s = line.strip().split('\t')
    protSpecies[p] = s
    
protPep = {}
f = open(peptideFile, 'rb')
test = 0
for line in f:
    if line.startswith('>'):
        if test == 1:
            protPep[name] = peptides
            peptides = []
            name = line.split(' ')[0][1:]
            name = name.split('.')[0]
        else:
            name = line.split(' ')[0][1:]
            name = name.split('.')[0]
            test = 1
            peptides = []
        
    else:
        if 'X' not in line:
            peptides.append(line.strip())
protPep[name] = peptides

speciesPeptides = {}
for key in protPep:
    species = '>' + protSpecies[key]
    peptides = protPep[key]
    speciesPeptides.setdefault(species, []).append(peptides)

for key in speciesPeptides:
    peptides = speciesPeptides[key]
    peptides = list(itertools.chain.from_iterable(peptides))
    speciesPeptides[key] = peptides


with open(output, 'wb') as out:
    for key in speciesPeptides:
        print >> out, key
        for pep in speciesPeptides[key]:
            print >> out, pep
