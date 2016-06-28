### this will take the .fa file and the dssp file and compare the two

import sys

fa= open(sys.argv[1]).read().splitlines()
dssp=open(sys.argv[2]).read().splitlines()

aminodic= {'A':1, 'R':2, 'N':3, 'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20} 
statedic={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2}
#0 is helix, 1 is strand, 2 is coil
amino= []


start=0



for line in dssp:
    value= line.split()
    if 'N-H-->O' in value:                       #indicates where to start parsing
        start= 1 
        continue                              #skip the header line
    if not start: 
        continue                              #skip all lines before the header line also
#### accessing the amino acids only #####
    if len(value[1])==4:
        amino.append(value[2])
        
    else:
        amino.append(value[3])

if len(amino) != len(fa[1]):
    print 'Warning!! Actual protein sequence != PDB prediction!!', len(fa[1]), len(amino)
else:
    print 'Actual protein sequence == PDB prediction'

count=0
for char, a in zip(fa[1], amino):
	count= count+1
	if char != a: 
		print 'missing or unmatching character: %s != %s' % char % a 
	if char== a:
		print char
