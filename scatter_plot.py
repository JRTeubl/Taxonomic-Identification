#!/usr/bin/env python
'''
This script visualizes the simulated, weighted results of 500 human and chimp samples. 
'''

import sys
import numpy as np
import pylab
import matplotlib.pyplot as plt
import collections
import itertools
from collections import Counter
import operator


seqFile1 = sys.argv[1]


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

f = open(seqFile1, 'rb')
speciesPeptides = openFile(f)


def weightPeptides(speciesPeptides):
    # all peptides are weighted such that the most unique peptides are given the most weight
    pepWeights = {}
    peptides = speciesPeptides.values()
    peptides = list(itertools.chain.from_iterable(peptides))
    totalPep = float(len(peptides))
    pepCounts = dict(Counter(peptides))
    pepCountsSorted = sorted(pepCounts.items(), key = operator.itemgetter(1))
    for x in pepCountsSorted:
        wp = float(1.0/x[1])/totalPep
        pepWeights[x[0]] = wp
    return pepWeights
pepWeights = weightPeptides(speciesPeptides)


def getPeptides(classificationsClass, speciesPeptides):
    classPeptides = {}
    for k in classificationsClass:
        if type(classificationsClass) is list:
            peptides = dividePeptides(k, speciesPeptides, speciesPeptides)
        else:
            peptides = dividePeptides(k, speciesPeptides, classificationsClass)
        if peptides:
            classPeptides[k] = sorted(peptides)
    return classPeptides

def dividePeptides(name, speciesPeptides, classifications):
    peptides = []
    try:
        sp = classifications[name]
        if speciesPeptides == classifications:
            for i in sp:
                peptides.append(i)
        else:
            for x in sp:
                try:
                    pep = speciesPeptides[x]
                    for i in pep:
                        peptides.append(i)
                except KeyError:
                    pass
    except KeyError:
        pass
    
    return peptides

def simulation(name, classPeptides):
    peptides = classPeptides[name]
    total = len(peptides)
    i = 0
    index = []
    while i < 500:
        index.append(np.random.randint(0,total, abs(np.random.normal(70, 40))))
        i = i + 1
    return index


def getSimPerSample(classPeptides):
    simSamples = {}
    for name in classPeptides:
        simPeptides = []
        index = simulation(name, classPeptides)
        pep = classPeptides[name]
        for i in index:
            temp = [pep[x] for x in i]
            temp = list(set(temp))
            simPeptides.append(temp)
        simSamples[name] = sorted(simPeptides)
    return simSamples


def commonPeptides(sim, pep, pepWeights):
    # a simulated data set is compared to a chosen species. The peptide weights of the shared peptides between both sets
    # are summed and the log of the sum is calculated for graphical use
    share = [x for x in sim if x in pep]
    shareWeight = np.log(np.sum([pepWeights[x] for x in share]) + .00001)
    return shareWeight

def calculateShare(sim, pep):
    cw = commonPeptides(sim, pep, pepWeights)
    return cw

classificationsMammals = {
    'primates' : ['Callithrix jacchus', 'Chlorocebus sabaeus', 'Gorilla gorilla', 'Homo sapiens', 'Macaca mulatta', 'Nomascus leucogenys',
                  'Pan troglodytes', 'Papio anubis', 'Tarsius syrichta'],
    'artiodactyla' : ['Bos taurus', 'Ovis aries', 'Sus scrofa', 'Tursiops truncatus', 'Vicugna pacos'],
    'carnivora' : ['Ailuropoda melanoleuca', 'Canis familiaris', 'Felis catus', 'Mustela putorius'],
    'chiroptera' : ['Myotis lucifugus', 'Pteropus vampyrus'],
    'eulipotyphla' : ['Erinaceus europaeus', 'Sorex araneus'],
    'lagomorpha' : ['Ochotona princeps', 'Oryctolagus cuniculus'],
    'rodentia' : ['Cavia porcellus', 'Dipodomys ordii', 'Ictidomys tridecemlineatus', 'Mus musculus', 'Rattus norvegicus'],
    'perissodactyla' : ['Equus caballus'],
    'afrosoricida' : ['Echinops telfairi'],
    'cingulata' : ['Dasypus novemcinctus'],
    'australidelphia' :['Sarcophilus harrisii', 'Monodelphis domestica', 'Macropus eugenii'],
    'monotremata' : ['Ornithorhynchus anatinus'],
    'scandentia': ['Tupaia belangeri']
}


primatePeptides = getPeptides(classificationsMammals['primates'], speciesPeptides)
humanPeps = primatePeptides['Homo sapiens']
gorillaPeps = primatePeptides['Gorilla gorilla']
macacaPeps = primatePeptides['Macaca mulatta']
papioPeps = primatePeptides['Papio anubis']
panPeps = primatePeptides['Pan troglodytes']
chlPeps = primatePeptides['Chlorocebus sabaeus']
calPeps = primatePeptides['Callithrix jacchus']
nomPeps = primatePeptides['Nomascus leucogenys']
tarPeps = primatePeptides['Tarsius syrichta']
    
simSamplesPrim = getSimPerSample(primatePeptides)
simHuman = simSamplesPrim['Homo sapiens']
hh = []
hg = []
hm = []
hpa = []
hpt = []
hcs = []
hcj = []
hnl = []
hts = []
for i in simHuman:
    x1 = calculateShare(i, humanPeps)
    hh.append(x1)
    x2 = calculateShare(i, gorillaPeps)
    hg.append(x2)
    x3 = calculateShare(i, macacaPeps)
    hm.append(x3)
    x4 = calculateShare(i, papioPeps)
    hpa.append(x4)
    x5 = calculateShare(i, panPeps)
    hpt.append(x5)
    x6 = calculateShare(i, chlPeps)
    hcs.append(x6)
    x7 = calculateShare(i, calPeps)
    hcj.append(x7)
    x8 = calculateShare(i, nomPeps)
    hnl.append(x8)
    x9 = calculateShare(i, tarPeps)
    hts.append(x9)

humanTest = {'Homo sapiens':hh, 'Gorilla gorilla':hg, 'Macaca mulatta':hm, 'Papio anubis':hpa, 'Pan troglodytes':hpt, 'Chlorocebus sabaeus': hcs,
             'Callithrix jacchus':hcj, 'Nomascus leucogenys': hnl, 'Tarsius syrichta':hts}

simPan = simSamplesPrim['Pan troglodytes']
ph = []
pg = []
pm = []
ppa = []
ppt = []
pcs = []
pcj = []
pnl = []
pts = []
for i in simPan:
    x1 = calculateShare(i, humanPeps)
    ph.append(x1)
    x2 = calculateShare(i, gorillaPeps)
    pg.append(x2)
    x3 = calculateShare(i, macacaPeps)
    pm.append(x3)
    x4 = calculateShare(i, papioPeps)
    ppa.append(x4)
    x5 = calculateShare(i, panPeps)
    ppt.append(x5)
    x6 = calculateShare(i, chlPeps)
    pcs.append(x6)
    x7 = calculateShare(i, calPeps)
    pcj.append(x7)
    x8 = calculateShare(i, nomPeps)
    pnl.append(x8)
    x9 = calculateShare(i, tarPeps)
    pts.append(x9)

panTest = {'Homo sapiens':ph, 'Gorilla gorilla':pg, 'Macaca mulatta':pm, 'Papio anubis':ppa, 'Pan troglodytes':ppt, 'Chlorocebus sabaeus': pcs,
             'Callithrix jacchus':pcj, 'Nomascus leucogenys': pnl, 'Tarsius syrichta':pts}


def scatterPlot(label, name):
    colors = ['Blue', 'BlueViolet', 'LightPink', 'Coral', 'Crimson', 'DarkGreen', 'DeepSkyBlue', 'Fuchsia', 'LimeGreen']
    c = 0
    xAxis = label[name]
    for i in label:
        yAxis = label[i]
        data = dict(zip(xAxis, yAxis))
        dataSorted = sorted(data.items(), key = lambda x:x[0], reverse = False)
        x = [k[0] for k in dataSorted]
        y = [k[1] for k in dataSorted]
        plt.scatter(x, y, color = colors[c], s = 10, alpha = .2)
        m, b = np.polyfit(x, y, 1)
        yNew = []
        for k in x:
            y1 = m*k+b
            yNew.append(y1)
        plt.plot(x, yNew, '-', c = colors[c], lw = 3, label = i)
        c = c + 1
    plt.legend(loc="upper left", fontsize = 8)
    plt.xlim(min(xAxis)-.01, max(xAxis)+.01)
    plt.ylim(min(xAxis)-.01, max(xAxis)+.01)
    plt.title('Feature Weights for %s Samples,\nsimulation n = 500'%name)
    plt.ylabel('%s log Fsj'%name, size = 12)
    plt.xlabel('log Fsj for all primates', size = 12)
    plt.savefig('/Users/jrteubl/Desktop/%s_scatter.pdf'%name)
    plt.close()

scatterPlot(humanTest, 'Homo sapiens')

scatterPlot(panTest, 'Pan troglodytes')
