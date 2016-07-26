from Bio.PDB import DSSP, PDBParser
import sys
from difflib import *
from sklearn.preprocessing import OneHotEncoder
import numpy as np

fa= open(sys.argv[1]).read().splitlines()
#dssp=open(sys.argv[2]).read().splitlines()    #this will eventually become the pdb file that we run dssp on
p=PDBParser()
structure = p.get_structure("1EFQA", "1EFQA.pdb")
model= structure[0]
dssp= DSSP(model, "1EFQA.pdb")
a_key = list(dssp.keys())#[31]
#print dssp[a_key]

for line in a_key:
	print dssp[line]




#### notes to self, for our structures there is only one model. So structure [0] is going to be the only one that works
#this is the new sliding table that I will add the the biopython dssp parser . Made 26 July








##############################################################################################################################
###################################################### SLIDING TABLE #########################################################
##############################################################################################################################

### get only the amino acids in file ## 
seq= fa[1]
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
