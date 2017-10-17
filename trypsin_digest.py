'''
Created on May 6, 2015

Input:
1. A fasta file with protein sequences on one line
2. The name and path for the output file. 
3. (Optional) a list of peptides with which to compared the in silico digested peptides. This can be used if there are only certain peptides
that should be considered, or if there are peptides that should be excluded (change line 51 to read: if x not in peptides)

Output:
1. The output file will contain all of the digested peptides for each protein. Each peptide will
be on its own line. Protein names will be the same as they were in the fasta file
'''

import sys



 
def trypsin(bases):
    sub = ''
    while bases:
        k, r = bases.find('K'), bases.find('R')
        cut = min(k, r)+1 if k > 0 and r > 0 else max(k, r)+1
        if cut == 0:
            sub += bases[:]
            bases = ''
        else:
            sub += bases[:cut]
            bases = bases[cut:]
        if not bases or bases[0] != 'P':
            yield sub
            sub = ''
def trypsinMissed(bases):
    sub = ''
    while bases:
        k,r = bases.find('K'), bases.find('R')
        cut1 = min(k,r)+1 if k>0 and r>0 else max(k,r)+1
        if cut1 == 0:
            sub += bases[:]
            bases = ''
        else:
            sub += bases[:cut1]
            bases = bases[cut1:]
            k,r = bases.find('K'), bases.find('R')
            cut2 = min(k,r)+1 if k>0 and r>0 else max(k,r)+1
            sub += bases[:cut2]
            bases = bases[cut2:]
        if not bases or bases[0] !='P':
            yield sub
            sub = ''

def main():

    filename = sys.argv[1]
    outname = sys.argv[2]
    if len(sys.argv) > 3:
        peptides = sys.argv[3]
    else:
        peptides = False
        
    out = open(outname, 'wb')
    if peptides:
        peps = open(peptides, 'rb')
        peptides = []
        for line in peps:
            peptides.append(line.split('\t')[0])
    
    f = open(filename, 'rb')
    for line in f:
        if line[0] == '>':
            out.write('%s'%line)
        else:
            p1 = list(trypsin(str(line.rstrip())))
            p2 = list(trypsinMissed(str(line.rstrip())))
            p = p1 + p2
            for x in p:
                if 40 >= len(x) >= 7:
                    if peptides:
                        if x in peptides:
                            out.write('%s\n'%x)
                    else:
                        out.write('%s\n'%x)

if __name__ == "__main__":
    main()    


