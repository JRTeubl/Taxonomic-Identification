#!/usr/bin/env python

'''
This script builds a linear SVM model using sequence date (seqFile1 and 2). It then test the model using experimental data imported as xml files.
MGF files are imported to confirm a chosen percentage of spectra are identified.

Classification are divided taxonomically: Class, order, species

Output: a table for each classification is made with normalized probabilities


*** make sure output includes version of python and scikit learn
import platform
print(platform.python_version())
print sklearn.__version__

'''

import numpy as np
import xml.etree.ElementTree as ET
import sys
import re
import os
import itertools
import string
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
style.use("ggplot")
from scipy import interp
import sklearn
import platform
from sklearn import svm, tree, neighbors
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from itertools import cycle
from collections import Counter
import operator


import platform
print(platform.python_version())
print sklearn.__version__

#Models:
#Class
#Order (mammal, fish, bird)
#Species (primate, carnivore, artidactyl, rodent, australidelphia)
#All models are an svm linear One vs rest classifier

seqFile1 = sys.argv[1]
seqFile2 = sys.argv[2]

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
speciesPeptides1 = openFile(f)
#speciesPeptides = openFile(f)

f = open(seqFile2, 'rb')
speciesPeptides2 = openFile(f)

# combine sequence files
speciesPeptides = {}
for key in speciesPeptides1:
    vals = speciesPeptides1[key]
    for v in vals:
        speciesPeptides.setdefault(key, []).append(v)
for key in speciesPeptides2:
    vals = speciesPeptides2[key]
    for v in vals:
        speciesPeptides.setdefault(key, []).append(v)

def weightPeptides(speciesPeptides):
    '''
    calculate peptide weights as the inverse of the number of times a peptide occurs divide by the total number of peptides
    '''
    pepWeights = {}
    peptides = speciesPeptides.values()
    peptides = list(itertools.chain.from_iterable(peptides))
    totalPep = float(len(peptides))
    pepCounts = dict(Counter(peptides))
    pepCountsSorted = sorted(pepCounts.items(), key = operator.itemgetter(1))
    for x in pepCountsSorted:
        wp1 = float(1.0/x[1])/totalPep
        #wp2 = float(1.0/x[1])/(totalPep - x[1])
        #wp3 = float(1.0/np.exp(x[1]))/totalPep
        pepWeights[x[0]] = wp1
    return pepWeights
pepWeights = weightPeptides(speciesPeptides)


def getPeptides(classificationsClass, speciesPeptides):
    '''
    for every classification obtain all peptides that map to a species within that classification
    '''
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
    '''
    extract peptides from each classification
    '''
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
    '''
    simulate 1000 sample from each species such that each sample contains a random list of peptides of varying length N(70,40)
    '''
    peptides = classPeptides[name]
    total = len(peptides)
    i = 0
    index = []
    while i < 1000:
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

def defineModel(simComp, label):
    simComp = np.asarray(simComp)
    data = simComp[:,0:-1]
    targets = simComp[:,-1]
    targets = label_binarize(targets, classes = label)
    X_train, X_test, y_train, y_test = train_test_split(data, targets, test_size=0.2)
    return X_train, y_train

svmMod = sklearn.multiclass.OneVsRestClassifier(svm.SVC(kernel='linear', probability=True))

classificationsClass = {
    'mammalia': ['Echinops telfairi', 'Bos taurus', 'Ovis aries', 'Sus scrofa', 'Tursiops truncatus', 'Vicugna pacos', 'Ailuropoda melanoleuca',
               'Canis familiaris', 'Felis catus', 'Mustela putorius', 'Myotis lucifugus', 'Pteropus vampyrus', 'Dasypus novemcinctus',
               'Sarcophilus harrisii', 'Monodelphis domestica', 'Macropus eugenii', 'Erinaceus europaeus', 'Sorex araneus','Ochotona princeps',
               'Oryctolagus cuniculus', 'Ornithorhynchus anatinus', 'Equus caballus', 'Callithrix jacchus','Chlorocebus sabaeus',
               'Gorilla gorilla', 'Homo sapiens', 'Macaca mulatta', 'Nomascus leucogenys', 'Pan troglodytes','Papio anubis',
               'Tarsius syrichta', 'Loxodonta africana', 'Cavia porcellus', 'Dipodomys ordii', 'Ictidomys tridecemlineatus','Mus musculus',
               'Rattus norvegicus', 'Tupaia belangeri', 'Lontra canadensis', 'Marmota monax'],
    'aves' : ['Anas platyrhynchos', 'Gallus gallus', 'Meleagris gallopavo', 'Ficedula albicollis', 'Taeniopygia guttata'],
    'actinoptergyii' : ['Oryzias latipes', 'Xiphophorus maculatus', 'Gadus morhua', 'Gasterosteus aculeatus', 'Lepisosteus oculatus',
                        'Oreochromis niloticus', 'Takifugu rubripes', 'Tetraodon nigroviridis'],
    'reptilia' : ['Anolis carolinensis', 'Pelodiscus sinensis']
} # current available species in ensmbl

classPeptides = getPeptides(classificationsClass, speciesPeptides)

if classPeptides:
    simSampleClass = getSimPerSample(classPeptides)
    simComp, labelClass = calculateShare(classPeptides, simSampleClass, pepWeights)
    X_trainClass, y_trainClass = defineModel(simComp, labelClass)


print 'finished class'

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
    'scandentia': ['Tupaia belangeri'],
    'proboscidea' : ['Loxodonta africana']
}

orderMamPeptides = getPeptides(classificationsMammals, speciesPeptides)

if orderMamPeptides:
    simSampleOrder = getSimPerSample(orderMamPeptides)
    simComp, labelMam = calculateShare(orderMamPeptides, simSampleOrder, pepWeights)
    X_trainMam, y_trainMam = defineModel(simComp, labelMam)

print "finished mammals"

classificationsAves = {
    'anseriforms' : ['Anas platyrhynchos'],
    'galliformes' : ['Gallus gallus', 'Meleagris gallopavo'],
    'passeriformes' : ['Ficedula albicollis', 'Taeniopygia guttata'],
}

orderAvesPeptides = getPeptides(classificationsAves, speciesPeptides)

if orderAvesPeptides:
    simSampleAvesOrder = getSimPerSample(orderAvesPeptides)
    simComp, labelAves = calculateShare(orderAvesPeptides, simSampleAvesOrder, pepWeights)
    X_trainAves, y_trainAves = defineModel(simComp, labelAves)

print "finished aves"

classificationsActin = {
    'beloniformes' : ['Oryzias latipes'],
    'cyprinodontiformes' : ['Xiphophorus maculatus'],
    'gadiforms' : ['Gadus morhua'],
    'gasterosteiformes' : ['Gasterosteus aculeatus'],
    'lepisosteiformes' : ['Lepisosteus oculatus'],
    'perciformes' : ['Oreochromis niloticus'],
    'tetraodontiformes' : ['Takifugu rubripes', 'Tetraodon nigroviridis']
}



orderActinPeptides = getPeptides(classificationsActin, speciesPeptides)
if orderActinPeptides:
    simSampleActinOrder = getSimPerSample(orderActinPeptides)
    simComp, labelActin = calculateShare(orderActinPeptides, simSampleActinOrder, pepWeights)
    X_trainActin, y_trainActin, = defineModel(simComp, labelActin)

print 'finished actinoptergyii'

primatePeptides = getPeptides(classificationsMammals['primates'], speciesPeptides)
if primatePeptides:
    simSamplesPrim = getSimPerSample(primatePeptides)
    simComp, labelPrim = calculateShare(primatePeptides, simSamplesPrim, pepWeights)
    X_trainPrim, y_trainPrim = defineModel(simComp, labelPrim)

print 'finished primates'

artPeptides = getPeptides(classificationsMammals['artiodactyla'], speciesPeptides)
if artPeptides:
    simSamplesArt = getSimPerSample(artPeptides)
    simComp, labelArt = calculateShare(artPeptides, simSamplesArt, pepWeights)
    X_trainArt, y_trainArt, = defineModel(simComp, labelArt)

print 'finished artiodactyla'


carnPeptides = getPeptides(classificationsMammals['carnivora'], speciesPeptides)
if carnPeptides:
    simSamplesCarn = getSimPerSample(carnPeptides)
    simComp, labelCarn = calculateShare(carnPeptides, simSamplesCarn, pepWeights)
    X_trainCarn, y_trainCarn = defineModel(simComp, labelCarn)

print 'finished carnivora'

rodentiaPeptides = getPeptides(classificationsMammals['rodentia'], speciesPeptides)
if rodentiaPeptides:
    simSamplesRod = getSimPerSample(rodentiaPeptides)
    simComp, labelRod = calculateShare(rodentiaPeptides, simSamplesRod, pepWeights)
    X_trainRod, y_trainRod = defineModel(simComp, labelRod)

print 'finished rodentia'


australiPeptides = getPeptides(classificationsMammals['australidelphia'], speciesPeptides)
if australiPeptides:
    simSamplesAus = getSimPerSample(australiPeptides)
    simComp, labelAus = calculateShare(australiPeptides, simSamplesAus, pepWeights)
    X_trainAus, y_trainAus = defineModel(simComp, labelAus)

print 'finished australidelphia'


def main():
    path = sys.argv[3]
    if len(sys.argv) > 4:
        cutoff = sys.argv[4]
        cutoff = float(cutoff)
    else:
        cutoff = .20        #cutoff for spectra coverage will default to 20% if no cutoff is included
    
    xmls = []
    mgfs = []
    gaps = []
    for filename in os.listdir(path):
        if filename.endswith('xml'):
            xmls.append(filename)
        elif filename.endswith('mgf'):
            mgfs.append(filename)
        else:
            print '%s is not mgf or xml file types'%filename
    
    classOut = pd.DataFrame(columns = ['actinoptergyii', 'aves', 'mammalia', 'reptilia'])
    
    orderMammalOut = pd.DataFrame(columns = ['artiodactyla', 'australidelphia', 'carnivora', 'chiroptera', 'cingulata',
                                            'eulipotyphla',  'lagomorpha', 'perissodactyla', 'primates', 'rodentia'])
   
    orderBirdOut = pd.DataFrame(columns = ['anseriforms', 'galliformes', 'passeriformes'])
    orderFishOut = pd.DataFrame(columns = ['beloniformes', 'cyprinodontiformes', 'gadiforms', 'gasterosteiformes', 'lepisosteiformes',
                                           'perciformes','tetraodontiformes'])
    speciesPrimateOut = pd.DataFrame(columns = ['Callithrix jacchus', 'Chlorocebus sabaeus', 'Gorilla gorilla', 'Homo sapiens','Papio anubis',
                                                'Nomascus leucogenys','Pan troglodytes', 'Macaca mulatta', 'Tarsius syrichta'])
    
    speciesArtOut = pd.DataFrame(columns = ['Bos taurus', 'Ovis aries', 'Sus scrofa', 'Tursiops truncatus', 'Vicugna pacos'])
    #speciesArtOut = pd.DataFrame(columns = ['Bos taurus', 'Sus scrofa', 'Tursiops truncatus', 'Vicugna pacos'])
    speciesCarnOut = pd.DataFrame(columns = ['Ailuropoda melanoleuca', 'Canis familiaris', 'Mustela putorius'])

    speciesRodOut = pd.DataFrame(columns = ['Cavia porcellus', 'Dipodomys ordii', 'Ictidomys tridecemlineatus', 'Mus musculus'])
    speciesAusOut = pd.DataFrame(columns = ['Sarcophilus harrisii', 'Monodelphis domestica', 'Macropus eugenii'])
    
    ### For col1a1 or col1a2 only ###
    # orderMammalOut = pd.DataFrame(columns = ['artiodactyla', 'australidelphia', 'carnivora', 'chiroptera', 'cingulata',
    #                                        'eulipotyphla',  'lagomorpha',  'monotremata', 'perissodactyla', 'primates', 'rodentia'])
    # orderMammalOut = pd.DataFrame(columns = ['australidelphia', 'carnivora', 'chiroptera',
    #                                       'eulipotyphla',  'lagomorpha',  'monotremata', 'perissodactyla', 'primates', 'rodentia'])

    # speciesPrimateOut = pd.DataFrame(columns = ['Callithrix jacchus', 'Chlorocebus sabaeus', 'Gorilla gorilla', 'Homo sapiens', 'Macaca mulatta',
    #                                          'Nomascus leucogenys'])
    # speciesPrimateOut = pd.DataFrame(columns = ['Callithrix jacchus', 'Homo sapiens', 'Pan troglodytes', 'Papio anubis', 'Tarsius syrichta'])
    # speciesCarnOut = pd.DataFrame(columns = ['Ailuropoda melanoleuca', 'Canis familiaris', 'Felis catus', 'Mustela putorius'])
    # speciesCarnOut = pd.DataFrame(columns = ['Ailuropoda melanoleuca', 'Felis catus'])
    # speciesRodOut = pd.DataFrame(columns = ['Cavia porcellus', 'Dipodomys ordii', 'Ictidomys tridecemlineatus', 'Mus musculus', 'Rattus norvegicus'])
    # speciesRodOut = pd.DataFrame(columns = ['Cavia porcellus', 'Mus musculus', 'Rattus norvegicus'])
    
    for xmlFile in xmls:
        name=re.sub('.xml','', xmlFile)
        for mgfFile in mgfs:
            mgfName = re.sub('.mgf', '', mgfFile)
            if name == mgfName:
                print name
                tree = ET.parse(path+xmlFile)
                root = tree.getroot()
                parentMap = {c:p for p in tree.iter() for c in p}
                                
                def removeWhitespaces(text):
                    if re.search('[a-zA-Z]', text):
                        sequence = text.translate(None, string.whitespace)
                    return sequence
                
                def peptideDict(root):
                    seqId = {}
                    for group in root.iter('group'):
                        if group.get('type') == 'model':
                            for dom in group.iter('domain'):
                                seq = dom.get('seq')
                        if group.get('type') == 'support':
                            for x in group.iter('note'):
                                if x.get('label') == 'Description':
                                    scanLine = x.text
                                    scanLine = removeWhitespaces(scanLine)
                                    scan = re.search('Locus:(.*)File', scanLine)
                                    scan = scan.group(1)
                                    seqId[(scan)] = (seq)
                    return seqId
                
                seqId = peptideDict(root)
            
                    
                def spectraDict(mgfFile):
                    spectraId = {}
                    f = open(path+mgfFile, 'rb')
                    for line in f:
                        if line.startswith('BEGIN IONS'):
                            mzs = []
                        if line.startswith('TITLE='):
                            scanLine = removeWhitespaces(line)
                            scan = re.search('Locus:(.*)File', scanLine)
                            scan = scan.group(1)
                        if line.startswith('PEPMASS='):
                            mh = float(line.strip('\r\n')[8:])
                        if line[0].isdigit():
                            mzs.append(line.split(' ')[0])
                        if line.startswith('END IONS'):
                            mzs = map(float, mzs)
                            spectraId[scan] = (mh, mzs)
                    return spectraId
                            
                spectraId = spectraDict(mgfFile)

                
                # First filtering step: check the number of identified spectra against the number of produced spectra. Discard
                # samples with very low ratios
                if float(len(seqId))/float(len(spectraId)) >= cutoff:
                    def peptideFragments(seqId, spectraId):
                        peptideFrag = {}
                        for k in seqId:
                            frags = spectraId[k]
                            peptideFrag[seqId[k]] = frags
                        return peptideFrag
                    
                    peptideFrags = peptideFragments(seqId, spectraId)
                    numPepB4filter =  len(peptideFrags)

                    
                    def getFalsePositives(root):
                        '''
                        If a decoy database is used, find all decoy hits and determine their evalues
                        '''
                        falseEvals = []
                        for group in root.iter('group'):
                            if group.get('type') == 'model':
                                evalue = group.get('expect')
                                for prot in group.iter('protein'):
                                    for x in prot.iter('note'):
                                        if x.get('label') == 'description':
                                            name = x.text
                                            if re.search(r'reversed', name):
                                                falseEvals.append(float(evalue))
                        return falseEvals
                    
                    falseEvals = getFalsePositives(root)

                    
                    def fdr(falseEvals, delta):
                        m = len(falseEvals)
                        evalsSort = sorted(falseEvals)
                        reject = []
                        for i in range(1, m):
                            test = np.multiply(np.divide(float(i),float(m)), delta)
                            if test < evalsSort[i]:
                                reject.append(evalsSort[i])
                        if reject:
                            threshold = np.min(reject)
                        elif falseEvals:
                            threshold = np.max(falseEvals)
                        else:
                            threshold = .01
                        return threshold
                    
                    threshold = fdr(falseEvals, .01)
                    
                    def removeBelowThreshold(root, threshold):
                        '''
                        add peptides with evalues below the threshold to the contaminants list
                        '''
                        pepEval = {}
                        contaminant = []
                        for group in root.iter('group'):
                            evalue = group.get('expect')
                            if evalue != None:
                                for domain in group.iter('domain'):
                                    seq = domain.get('seq')
                                    pepEval.setdefault(seq, []).append(evalue)
                            for k in pepEval:
                                val = min(map(float,pepEval[k]))
                                if val >= threshold:
                                    contaminant.append(k)
                            contaminant = list(set(contaminant))
                        return contaminant
                    
                    contaminant = removeBelowThreshold(root, threshold)
                    
                    def removeFromCrap(root, parentMap, contaminant):
                        '''
                        create a list of peptides that map to proteins in the crapome
                        '''
                        for fileXml in root.iter('file'):
                            db = fileXml.get('URL')
                            if db.endswith('crap.fasta'):
                                parentProtein = parentMap[parentMap[fileXml]]
                                for domain in parentProtein.iter('domain'):
                                    seq = domain.get('seq')
                                contaminant.append(seq)
                            contaminant = list(set(contaminant))
                        return contaminant
                    
                    contaminant = removeFromCrap(root, parentMap, contaminant)
                    #print len(contaminant)
                    
                    def removeHydPro(root, parentMap, contaminant):
                        '''
                        Hydroxy prolines can only exist at the Y residue of the pattern G-X-Y.
                        This function removes peptides with hydroxy prolines in any other position
                        '''
                        cont = []
                        tempKeep = []
                        for domain in root.iter('domain'):
                            seq = domain.get('seq')
                            start = domain.get('start')
                            if seq.startswith('G'):
                                aas = []
                                for aa in domain.iter('aa'):
                                    aminoacid = aa.get('type')
                                    if aminoacid == 'P':
                                        at = aa.get('at')
                                        placement = float(at) - float(start) + 1
                                        aas.append(placement)
                                tempCont = []
                                for x in aas:
                                    if x%3.0 != 0:
                                        tempCont.append(seq)
                                if len(tempCont) == 0:
                                    tempKeep.append(seq)
                                else:
                                    cont.append(seq)
                        for x in cont:
                            if x not in tempKeep:
                                contaminant.append(x)
                        return contaminant
        
                    contaminant = removeHydPro(root, parentMap, contaminant)
                    
                    def removeHighMissedCleavage(root, parentMap, contaminant):
                        '''
                        find spectra matches with more than one missed cleavage and remove them. 
                        '''
                        for domain in root.iter('domain'):
                            numMissClev = int(domain.get('missed_cleavages'))
                            if numMissClev > 1:
                                seq = domain.get('seq')
                                contaminant.append(seq)
                        return contaminant
                    
                    contaminant = removeHighMissedCleavage(root, parentMap, contaminant)
                    contaminant = list(set(contaminant))

                    
                    def removeContaminants(contaminant, peptideFrags):
                        for key in contaminant:
                            if key in peptideFrags:
                                del peptideFrags[key]
                        return peptideFrags
                    
                    cleanPepFrag = removeContaminants(contaminant, peptideFrags)
                    peptides = cleanPepFrag.keys()
                    numPepAfterfilter =  len(peptides)
                    
                    outputPath = path+'filters/MC20perCutoff/'
                    outputName = outputPath + name+'-'+str(cutoff)+'-'+str(numPepB4filter)+'-'+str(numPepAfterfilter)
                    if not os.path.exists(outputPath):
                        os.mkdir(outputPath)
                    if not os.path.exists(outputName):
                        os.mkdir(outputName)
                        
                    
                    def getPrb(levelPeptides, peptides, x_train, y_train):
                        sampleClass, label = calculateShare(levelPeptides, peptides, pepWeights)
                        sampleClass = np.asarray(sampleClass)
                        sampleClass = np.reshape(sampleClass, (1,np.product(sampleClass.shape)))
                        ClassMod =  svmMod.fit(x_train, y_train)
                        prob = ClassMod.predict_proba(sampleClass)
                        prob = list(prob[0])
                        return prob, label
                    
                    def outputAndSelect(levelPeptides, x_train, y_train, level, label):
                        prob, x = getPrb(levelPeptides, peptides, x_train, y_train)
                        data = dict(zip(label, prob))
                        probdf = pd.DataFrame(data, index = [name])
                        probdf.to_csv(outputName+'/'+name+'_'+level+'.txt', sep='\t')
                        ix = prob.index(max(prob))
                        select = label[ix]
                        # if max(prob) >= .75:
                        #     select = label[ix]
                        # else:
                        #     select = None
                        return select, probdf

                    selectClass, classDF = outputAndSelect(classPeptides, X_trainClass, y_trainClass, 'class', labelClass)
                    #print selectClass
                    classOut = classOut.append(classDF)
                    if selectClass == 'mammalia':
                        selectOrder, orderDF = outputAndSelect(orderMamPeptides, X_trainMam, y_trainMam, 'order', labelMam)
                        #print selectOrder
                        orderMammalOut= orderMammalOut.append(orderDF)
                        if selectOrder == 'primates':
                            selectPrim, speciesDF = outputAndSelect(primatePeptides, X_trainPrim, y_trainPrim, 'primate', labelPrim)
                            #print selectPrim
                            speciesPrimateOut = speciesPrimateOut.append(speciesDF)
                        elif selectOrder == 'carnivora':
                            selectCarn, speciesDF = outputAndSelect(carnPeptides, X_trainCarn, y_trainCarn, 'carnivora', labelCarn)
                            #print selectCarn
                            speciesCarnOut = speciesCarnOut.append(speciesDF)
                        elif selectOrder == 'rodentia':
                            selectRod, speciesDF = outputAndSelect(rodentiaPeptides, X_trainRod, y_trainRod, 'rodentia', labelRod)
                            #print selectRod
                            speciesRodOut = speciesRodOut.append(speciesDF)
                        elif selectOrder == 'artiodactyla':
                           selectArt, speciesDF = outputAndSelect(artPeptides, X_trainArt, y_trainArt, 'artidactyla', labelArt)
                           #print selectArt
                           speciesArtOut = speciesArtOut.append(speciesDF)
                        elif selectOrder == 'australidelphia':
                            selectAus, speciesDF = outputAndSelect(australiPeptides, X_trainAus, y_trainAus, 'australidelphia', labelAus)
                            #print selectAus
                            speciesAusOut = speciesAusOut.append(speciesDF)
                        elif selectOrder == None:
                            print 'The order of sample %s could not be determined with greater than 75 percent certainty.'%name
                        
                    if selectClass == 'aves':
                        selectAves, classDF = outputAndSelect(orderAvesPeptides, X_trainAves, y_trainAves, 'order', labelAves)
                        #print selectAves
                        orderBirdOut = orderBirdOut.append(classDF)
                    if selectClass == 'actinoptergyii':
                        selectAct, classDF = outputAndSelect(orderActinPeptides, X_trainActin, y_trainActin, 'order', labelActin)
                        #print selectAct
                        orderFishOut = orderFishOut.append(classDF)
                    if selectClass == None:
                        print 'The Class of sample %s could not be determined with greater than 75 percent certainty'%name
                

                else:
                    print 'sample has less than %.2f spectra coverage'%cutoff
                
                def normalizeToOne(df):
                    df['sum'] = df.sum(axis = 1)
                    sums = df.pop('sum')
                    df = df.div(sums, axis = 0)
                    return df
                
                classOut = normalizeToOne(classOut)
                classOut.to_csv(outputPath+'class.txt', sep = '\t')
                
                orderBirdOut = normalizeToOne(orderBirdOut)
                orderBirdOut.to_csv(outputPath+'orderBird.txt', sep = '\t')
                orderMammalOut = normalizeToOne(orderMammalOut)
                orderMammalOut.to_csv(outputPath+'orderMammal.txt', sep = '\t')
                orderFishOut = normalizeToOne(orderFishOut)
                orderFishOut.to_csv(outputPath+'orderFish.txt', sep = '\t')
                
                speciesPrimateOut = normalizeToOne(speciesPrimateOut)
                speciesPrimateOut.to_csv(outputPath+'speciesPrimate.txt', sep = '\t')
                
                speciesArtOut = normalizeToOne(speciesArtOut)
                speciesArtOut.to_csv(outputPath+'speciesArt.txt', sep = '\t')
                speciesCarnOut = normalizeToOne(speciesCarnOut)
                speciesCarnOut.to_csv(outputPath+'speciesCarn.txt', sep = '\t')
                speciesRodOut = normalizeToOne(speciesRodOut)
                speciesRodOut.to_csv(outputPath+'speciesRod.txt', sep = '\t')
                speciesAusOut = normalizeToOne(speciesAusOut)
                speciesAusOut.to_csv(outputPath+'speciesAus.txt', sep = '\t')
                
if __name__ == "__main__":
    main()

