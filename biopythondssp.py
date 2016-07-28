from Bio.PDB import DSSP, PDBParser
import sys
from difflib import *
from sklearn.preprocessing import OneHotEncoder
import numpy as np


##############################################################################################################################
###################################################### PARSE DSSP ############################################################
##############################################################################################################################

fa= open(sys.argv[1]).read().splitlines()
filename=sys.argv[2]    #this will eventually become the pdb file that we run dssp on
dssp=open(filename)
filename= filename[-9:]     #only the actual filename t.ex 1EFQA.pdb
p=PDBParser()
structure = p.get_structure(filename, dssp)
model= structure[0]        #there is only one structure for dssp (NMR for example has more) and the dssp parser can only take one structure
dssp= DSSP(model, filename)
a_key = list(dssp.keys())

statedic={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2, '-':'XXX'}
#0 is helix, 1 is strand, 2 is coil
dsspAA=[ ]
states= [ ] 
for line in a_key:
	print dssp[line]
	dsspAA.append(dssp[line][1])
	states.append(statedic[dssp[line][2]])
print states
##############################################################################################################################
######################################################ONE HOT ENCODING #######################################################
##############################################################################################################################

seq= fa[1]
#dsspAA seen above in parser
d=Differ()
diff= d.compare(seq, dsspAA)
comp= '\n'.join(diff)


#states=np.asarray(states).reshape(-1,1)
NumsOnly= []
final= [ ] 
for val in states:
	if val == 'XXX':
		print '0 , 0, 0'
		final.append(' 0 0 0 ')
	else:
		NumsOnly.append(val)
		NumsOnlyArray=np.asarray(NumsOnly).reshape(-1,1)
		enc= OneHotEncoder()
		encoded= np.around(enc.fit_transform(NumsOnlyArray).toarray(), decimals= 0)
		encoded= encoded.tolist()
print encoded





##############################################################################################################################
###################################################### SLIDING TABLE #########################################################
##############################################################################################################################

#seq= fa[1] (just a reminder what seq is- appears in OHE)
zeros= 15        #change this accordingly to how many zeros you want on each side. Be sure to think of how many lines of zeros you will have in the OHE part
seq= ('o' * zeros) + seq + ('o' * zeros)
print seq

### iterate through each character in seq and make a row per character. Each row has 1 one and the rest are zeros 
letters= 'ACDEFGHIKLMNPQRSTVWY'   # to create numbers for each amino acid
def encode (seq, l):
    row= []
    for c in seq:
        if c == l:
            row.append(1)
        else:
            row.append(0)
    return row

row= []
for l in letters:
    row.append(encode(seq, l))
WholeSeq= np.asarray(row)

final=[]
for index in xrange(len(seq)):
		slide= WholeSeq[0:,index:index+zeros] 				#0: == use all rows ',' <column no> ':' <column no>
		if len(slide[0]) == zeros:
			pass
#			print slide
			for line in slide: 					# because you cant append arrays to each other, change to a list and then back to array. Doesn't seem to hurt the values at all
				final.append(line)
final=np.asarray(final)
#print final
out=open('blabblabh','w')
np.savetxt(out, np.around(final, decimals=0), fmt='%.0f')
out.close()
print 'Sliding table found at: ', out
