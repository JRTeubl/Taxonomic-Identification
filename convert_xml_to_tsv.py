#!/usr/bin/env python
'''
This script converts an xml file to a tsv file with the following headers: file, scan, charge, sequence.
The sequence includes the modifications in the format: ELAEDGC[+57.0]SGVEVR

'''

import os
import sys
import xml.etree.ElementTree as ET
import re
from string import whitespace
import pandas as pd
import csv

path = sys.argv[1]
output = sys.argv[2]

for filename in os.listdir(path):
    if filename.endswith('.xml'):
        name = filename.split('.')[0]
        tree = ET.parse(path+filename)
        root = tree.getroot()

        
        def removeWhitespaces(text):
            if re.search('[a-zA-Z]', text):
                sequence = text.translate(None, whitespace)
            return sequence
    
        def groupDict(root):
            groupDict = {}
            for group in root.iter('group'):
                if group.get('type') == 'model':
                    iden = group.get('id')
                    charge = group.get('z') 
                if group.get('type') == 'support':
                    for x in group.iter('note'):
                        if x.get('label') == 'Description':
                            scanLine = x.text
                            scanLine = removeWhitespaces(scanLine)
                            scan = re.search('Locus:(.*)File', scanLine)
                            scan = scan.group(1)
                            groupDict[iden] = [charge, scan]
            return groupDict
        
        groupDict = groupDict(root)
        
        def addMod(place, seq):
            l = len(place)
            p = [x[0] for x in place]
            m = [x[1] for x in place]
            if l > 0:
                while l > 0:
                    seq = seq[0:p[l-1]]+m[l-1]+seq[p[l-1]:]
                    l = l-1

            seq = str(seq)
            return seq
        
        def domainDict(root):
            domainDict = {}
            for domain in root.iter('domain'):
                modified = []
                iden = domain.get('id')
                idsub = iden.split('.')[0]
                seq = domain.get('seq')
                start = domain.get('start')
                place = []
                for aa in domain.iter('aa'):
                    at = aa.get('at')
                    mass = aa.get('modified')
                    # if not mass.startswith('-'):
                    #     mass = '(+' + mass + ')'
                    # else:
                    #     mass = '(' + mass + ')'
                    placement = int(float(at) - float(start) + 1)
                    place.append([placement, str(mass)])
                place = sorted(place)
                #seqMod = addMod(place, seq)
                
                domainDict[iden] = seq
            return domainDict
        
        domainDict = domainDict(root)

        
        
        forDf = []

        for iden in domainDict:
            seq= domainDict[iden]
            idenSplit = iden.split('.')
            protId = '.'.join((idenSplit[0],idenSplit[1]))
            groupId = idenSplit[0]
            charge, scan = groupDict[groupId]

            forDf.append([name, scan, charge, seq])
            
        dfPeptides = pd.DataFrame(forDf, columns = ['file','scan', 'charge','sequence'])
        
        with open(output+name+'.tsv', 'wb') as out:
            dfPeptides.to_csv(out, sep='\t', quoting=csv.QUOTE_NONE, index=False) 
            
            
            
            
            
