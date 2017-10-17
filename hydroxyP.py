#!/usr/bin/env python
 
 '''
 This script identifies peptides with incorrect hydroxylated prolines. Hydroylated prolines can only be in the
 Y position of the GXY pattern found in collagen.
 '''
import sys
import re
import xml.etree.ElementTree as ET

xmlFile = sys.argv[1]
name=re.sub('.xml','', xmlFile)
tree = ET.parse(xmlFile)
root = tree.getroot()
parentMap = {c:p for p in tree.iter() for c in p}

contaminant = []

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

removeHydPro(root, parentMap, contaminant)
print contaminant
