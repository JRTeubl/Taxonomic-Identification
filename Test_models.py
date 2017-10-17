#!/usr/bin/env python
'''

This script test 4 different machine learning algorithms: SVM linear, SVM polynomial, Nearest neighbor, and decision tree. It outputs an AUC curve
comparing the average score for classification for each model. All parsing functions are described in Taxonomic_Model.py
'''


import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
from sklearn import svm, tree, neighbors
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from scipy import interp
import pylab
import itertools
from collections import Counter
import operator
import time


seqFile1 = sys.argv[1]
#seqFile2 = sys.argv[2]
output = sys.argv[2]

classificationsMammals = {
    'primates' : ['Callithrix jacchus', 'Chlorocebus sabaeus', 'Gorilla gorilla', 'Homo sapiens', 'Macaca mulatta', 'Nomascus leucogenys',
                  'Pan troglodytes', 'Papio anubis', 'Tarsius syrichta'],
    'artiodactyla' : ['Bos taurus', 'Ovis aries', 'Sus scrofa', 'Tursiops truncatus', 'Vicugna pacos', 'Pantholops hodgsonii'],
    'carnivora' : ['Ailuropoda melanoleuca', 'Canis familiaris', 'Felis catus', 'Mustela putorius'],
    'chiroptera' : ['Myotis lucifugus', 'Pteropus vampyrus'],
    'eulipotyphla' : ['Erinaceus europaeus', 'Sorex araneus'],
    'lagomorpha' : ['Ochotona princeps', 'Oryctolagus cuniculus'],
    'rodentia' : ['Cavia porcellus', 'Dipodomys ordii', 'Ictidomys tridecemlineatus', 'Mus musculus', 'Rattus norvegicus'],
    'perissodactyla' : ['Equus caballus'],
    'afrosoricida' : ['Echinops telfairi'],
    'cingulata' : ['Dasypus novemcinctus'],
    'dasyuromorphia' : ['Sarcophilus harrisii'],
    'didelphimophia' : ['Monodelphis domestica'],
    'diprotodontia' : ['Macropus eugenii'],
    'monotremata' : ['Ornithorhynchus anatinus'],
    'scandentia': ['Tupaia belangeri']
}

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
#speciesPeptides1 = openFile(f)
speciesPeptides = openFile(f)

#f = open(seqFile2, 'rb')
#speciesPeptides2 = openFile(f)

# speciesPeptides = {}
# for key in speciesPeptides1:
#     vals = speciesPeptides1[key]
#     for v in vals:
#         speciesPeptides.setdefault(key, []).append(v)
# for key in speciesPeptides2:
#     vals = speciesPeptides2[key]
#     for v in vals:
#         speciesPeptides.setdefault(key, []).append(v)

def weightPeptides(speciesPeptides):
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
    while i < 1000:
        index.append(np.random.randint(0,total, abs(np.random.normal(70,40))))
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
    share = [x for x in sim if x in pep]
    shareWeight = np.sum([pepWeights[x] for x in share])
    return shareWeight

def calculateShare(levelPeptides, simSampleLevel, pepWeights):
    simComp = []
    label = []
    if type(simSampleLevel) is dict:
        for k in simSampleLevel:
            spSim = simSampleLevel[k]
            for sim in spSim:
                sharedPep = []
                for i in levelPeptides:
                    pep = levelPeptides[i]
                    pw = np.sum([pepWeights[x] for x in pep])
                    cw = commonPeptides(sim, pep, pepWeights)
                    share = cw/pw
                    #print '%f, %f'%(cw,pw)
                    sharedPep.append(share)   
                sharedPep.append(k)
                simComp.append(sharedPep)
            label.append(k)
    else:
        for i in levelPeptides:
            pep = levelPeptides[i]
            pw = np.sum([pepWeights[x] for x in pep])
            cw = commonPeptides(simSampleLevel, pep, pepWeights)
            share = cw/pw
            simComp.append(share)
            label.append(i)
    return simComp, label

primatePeptides = getPeptides(classificationsMammals['primates'], speciesPeptides)


simSamplesPrim = getSimPerSample(primatePeptides)
simComp, label = calculateShare(primatePeptides, simSamplesPrim, pepWeights)
simComp = np.asarray(simComp)
data = simComp[:,0:-1]
targets = simComp[:,-1]
targets = label_binarize(targets, classes = label)
n_classes = targets.shape[1]

#Models to test
svmLin = OneVsRestClassifier(svm.SVC(kernel='linear', probability=True))
svmPoly = OneVsRestClassifier(svm.SVC(kernel='poly', probability=True))
dtc = tree.DecisionTreeClassifier()
k = 15
nnk = neighbors.KNeighborsClassifier(k)

def macroROC(fpr, tpr, n_classes):
    '''
    computes the average score for all classes in a classification method
    '''
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(n_classes)]))
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(n_classes):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    mean_tpr /= n_classes
    mean_auc = auc(all_fpr, mean_tpr)
    return mean_tpr, all_fpr, mean_auc

def ROC(data, targets, classifier, n_classes):
    '''
    calculates true positives and false positives and area under the curve 
    '''
    fpr = {}
    tpr = {}
    roc_auc= {}
    X_train, X_test, y_train, y_test = train_test_split(data, targets, test_size=0.75)
    if hasattr(classifier, 'decision_function'):
        y_score = classifier.fit(X_train, y_train).decision_function(X_test)
    else:
        y_score = classifier.fit(X_train, y_train).predict(X_test)
    for i in range(n_classes):
        fpr[i], tpr[i], threshold = roc_curve(y_test[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    return tpr, fpr, roc_auc

def rocPlot(fprs, tprs, aucs, names, colors, title, saveAs):
    '''
    plot roc curve
    '''
    plt.figure()
    for i in range(len(names)):
        plt.plot(fprs[i], tprs[i], label='%s (area = %0.2f)'%(names[i], aucs[i]), color=colors[i], linestyle='--', linewidth=3)
    
    plt.plot([0, 1], [0, 1], 'k--', lw=2)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right", fontsize = 10)
    plt.savefig(output+saveAs+'.pdf')

start = time.time()
SVMlinTPR, SVMlinFPR, SVMlinAUC = ROC(data, targets, svmLin, n_classes)
end = time.time()
print(end - start)

start = time.time()
SVMpolyTPR, SVMpolyFPR, SVMpolyAUC = ROC(data, targets, svmPoly, n_classes)
end = time.time()
print(end - start)

start = time.time()
dtcTPR, dtcFPR, dtcAUC = ROC(data, targets, dtc, n_classes)
end = time.time()
print(end - start)

start = time.time()
nnkTPR, nnkFPR, nnkAUC = ROC(data, targets, nnk, n_classes)
end = time.time()
print(end - start)

cm = pylab.get_cmap('gist_rainbow')
colors = []
for i in range(len(label)):
    colors.append(cm(1.*i/len(label)))
rocPlot(SVMlinFPR, SVMlinTPR, SVMlinAUC, label, colors,'Support Vector Classifier - Linear', 'svm_lin')
rocPlot(SVMpolyFPR, SVMpolyTPR, SVMpolyAUC, label, colors, 'Support Vector Classifier - Poly', 'svm_poly')
rocPlot(dtcFPR, dtcTPR, dtcAUC, label, colors, 'Decision Tree Classifier', 'dtc')
rocPlot(nnkFPR, nnkTPR, nnkAUC, label, colors, 'Nearest Neighbor k=%i'%k, 'nnk')

SVMlinTPRavg, SVMlinFPRavg, SVMlinAUCavg = macroROC(SVMlinFPR, SVMlinTPR, n_classes)
SVMpolyTPRavg, SVMpolyFPRavg, SVMpolyAUCavg = macroROC(SVMpolyFPR, SVMpolyTPR, n_classes)
dtcTPRavg, dtcFPRavg, dtcAUCavg = macroROC(dtcFPR, dtcTPR, n_classes)
nnkTPRavg, nnkFPRavg, nnkAUCavg = macroROC(nnkFPR, nnkTPR, n_classes)
fprs = [SVMlinFPRavg, SVMpolyFPRavg, dtcFPRavg, nnkFPRavg]
tprs = [SVMlinTPRavg, SVMpolyTPRavg, dtcTPRavg, nnkTPRavg]
aucs = [SVMlinAUCavg, SVMpolyAUCavg, dtcAUCavg, nnkAUCavg]

models = ['SVM linear', 'SVM poly', 'Decision Tree', 'Nearest neighbor']
colorsAvgs = ['CornflowerBlue', 'Crimson', 'DarkOrange', 'DarkOliveGreen']
rocPlot(fprs,tprs,aucs, models, colorsAvgs, 'Average Multi-class Classifiers', 'average_classifiers')
