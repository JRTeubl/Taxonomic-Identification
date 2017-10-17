#!/usr/bin/env python

'''
This script outputs a boxplot describing the results of the svm analysis from Taxonomic_Model.py. Results are split by correct and incorrect
responses. For each sample there is one correct response and multiple incorrect responses. Each taxonomic class has it's own bar plot comparing
correct and incorrect probabilities.

'''

import sys
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

path = sys.argv[1]
classFile = sys.argv[2]
orderFile = sys.argv[3]
speciesFile = sys.argv[4]


def openAnswerKey(answerKey):
    f = open(answerKey)
    key = {}
    for line in f:
        n, c = line.strip('\n').split('\t')
        key[n] = c
    return key

classKey = openAnswerKey(classFile)
orderKey = openAnswerKey(orderFile)
speciesKey = openAnswerKey(speciesFile)

for filename in os.listdir(path):
    if filename == 'class.txt':
        dfClass = pd.read_table(path+filename, index_col = 0)
    if filename == 'orderMammal.txt':
        dfOrderMammal = pd.read_table(path+filename, index_col = 0)
    if filename == 'speciesArt.txt':
        dfSpeciesArt = pd.read_table(path+filename, index_col = 0)
    if filename == 'speciesCarn.txt':
        dfSpeciesCarn = pd.read_table(path+filename, index_col = 0)
    if filename == 'speciesPrimate.txt':
        dfSpeciesPrimate = pd.read_table(path+filename, index_col = 0)
    if filename == 'speciesRod.txt':
        dfSpeciesRod = pd.read_table(path+filename, index_col = 0)

def dataForPlot(df, answerKey):
    correct = []
    incorrect = []
    for i in df.index:
        temp = []
        try:
            cor = answerKey[i]
            if cor in df.columns:
                c = df.ix[i][cor]
                correct.append(c)
            else:
                correct.append(0.0)
            inc = []
            for col in df.columns:
                if col != cor:
                    nc = df.ix[i][col]
                    inc.append(nc)
            #incorrect.append(np.mean(inc))
            temp.append(inc)
        except KeyError:
            print i
        for x in temp:
            for i in x:
                incorrect.append(i)
    return correct, incorrect

classCor, classIncor = dataForPlot(dfClass, classKey)
mammalCor, mammalIncor = dataForPlot(dfOrderMammal, orderKey)
artCor, artIncor = dataForPlot(dfSpeciesArt, speciesKey)
carnCor, carnIncor = dataForPlot(dfSpeciesCarn, speciesKey)
primCor, primIncor = dataForPlot(dfSpeciesPrimate, speciesKey)
rodCor, rodIncor = dataForPlot(dfSpeciesRod, speciesKey)

colors = ['Crimson', 'CornflowerBlue', 'Coral', 'BlueViolet', 'LightGreen', 'LightSeaGreen']
labels = ['Class', 'Mammal', 'Artiodactyla', 'Carnivore', 'Primate', 'Rodent']

fig, axes = plt.subplots(nrows=1, ncols=6, figsize=(10, 6), sharey=True)

def fillBox(box, i, color, ax, data):
        bx = box['boxes'][i]
        boxX = []
        boxY = []
        for j in range(5):
            boxX.append(bx.get_xdata()[j])
            boxY.append(bx.get_ydata()[j])
        boxCoords = list(zip(boxX, boxY))
        boxPolygon = Polygon(boxCoords, facecolor=color)
        ax.add_patch(boxPolygon)
        med = box['medians'][i]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(med.get_xdata()[j])
            medianY.append(med.get_ydata()[j])
            ax.plot(medianX, medianY, 'k', c = 'black')
        ax.plot([np.average(med.get_xdata())], [np.average(data)],
                 color='black', marker='*', markeredgecolor='black')
        print np.average(data)
        
def plotBox(cor, incor, ax, color):
    box = ax.boxplot([cor, incor], widths = .7)
    plt.setp(box['boxes'], color='black')
    plt.setp(box['whiskers'], color='black')    
    fillBox(box, 0, color, ax, cor)
    fillBox(box, 1, color, ax, incor)    
    

plotBox(classCor, classIncor, axes[0], colors[0])
plotBox(mammalCor, mammalIncor, axes[1], colors[1])
plotBox(artCor, artIncor, axes[2], colors[2])
plotBox(carnCor, carnIncor, axes[3], colors[3])
plotBox(primCor, primIncor, axes[4], colors[4])
plotBox(rodCor, rodIncor, axes[5], colors[5])

i = 0
while i < len(labels):
    axes[i].set_xticklabels(['Correct', 'Incorrect']*6, rotation = 50, fontsize = 10)
    axes[i].set_title(labels[i])
    i = i + 1

plt.figtext(0.929, 0.711, ' Average\n Value', color='black', weight='roman',
            size='small', backgroundcolor='white')
plt.figtext(0.912, 0.719, '*', color='black', backgroundcolor='white',
            weight='roman', size='x-large')

plt.savefig(path+'final_fig.pdf')
#plt.show()
    


