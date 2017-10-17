#!/usr/bin/env python

'''
This script makes mgf files useable in Lorikeet by removing the last column of data
'''

import os
import re
import sys
folder = sys.argv[1] # folder where mgf files are
newFolder = sys.argv[2] # make new folder for new mgf files
if not os.path.exists(newFolder):
    newFolder = os.mkdirs(newFolder)
for f in os.listdir(folder):
    if f.endswith('.mgf'):
        fout = open(newFolder+f, 'wb')
        fop = open(folder+f, 'rb')
        for line in fop:
            if line[0].isdigit():
                line = line.strip().split(' ')[0:2]
                print >> fout, '%f %f'%(float(line[0]), float(line[1]))
            else:
                print >> fout, line.strip()
            