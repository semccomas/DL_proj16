from Bio.PDB import DSSP, PDBParser 
import sys 
from difflib import * 
from sklearn.preprocessing import OneHotEncoder, Imputer 
import numpy as np


##############################################################################################################################
###################################################### PARSE DSSP ############################################################
##############################################################################################################################

fa= open(sys.argv[1]).read().splitlines()
filename=sys.argv[2]    #this will eventually become the pdb file that we run dssp on
dssp=open(filename)
#### OBS: sys.argv[3] is OHE table, sys.argv[4] is sliding table

p=PDBParser()
structure = p.get_structure(filename, dssp)
model= structure[0]        #there is only one structure for dssp (NMR for example has more) and the dssp parser can only take one structure
dssp= DSSP(model, filename)
a_key = list(dssp.keys())

#statedic={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2, '-':np.nan}
#0 is helix, 1 is strand, 2 is coil
#above leads to only 3 OHE values when using imputer, below leads to 4. Leaving it like that because it seems to make more sense that it should have its own value if missing 
statedic={'H':60, 'I':60 , 'G':60, 'E':80, 'B':80, 'T':1, 'S':1, 'L':1, '-':np.nan}

dsspAA=[ ]
states= [ ] 
for line in a_key:
#	print dssp[line]
	dsspAA.append(dssp[line][1])
	states.append(statedic[dssp[line][2]])

##############################################################################################################################
######################################################ONE HOT ENCODING #######################################################
##############################################################################################################################

seq= fa[1]
#dsspAA seen above in parser
d=Differ()
diff= d.compare(seq, dsspAA)
comp= '\n'.join(diff)
print states, len(states)
print comp, comp.count('\n')
match=0
total= [ ]
for diff in comp.split('\n'):
	if '-' not in diff and 'X' not in diff:
		total.append(states[match])
		match= match+1
	elif 'X' in diff:
		match= match + 1
	else:
		total.append(np.nan)

#### adding zeros to pad for the padded slider:
n= np.nan
Fpad= [n, n, n, n, n, n, n, n]      #using 9 values here because you want to start OHE for real when the feature is in the middle of the table
padded = Fpad + total + Fpad
#print padded, len(padded)

###


total=np.asarray(padded).reshape(-1,1)
imputer = Imputer(missing_values='NaN', strategy='mean').fit(total)
imputer = imputer.transform(total)

enc=OneHotEncoder()
encoded=np.around(enc.fit_transform(imputer).toarray(), decimals=0)

minus= len(seq)-len(dsspAA)

OHE=open(sys.argv[3], 'w')
np.savetxt(OHE, np.around(encoded, decimals=0), fmt='%.0f')
#print encoded, len(encoded)
#print total, len(total)
#print states, 'states'
print


print 'The length of this file is now: %s . Should == the length of the sequence+ 16 -->' % len(encoded) ,len(seq), 'which is:', (len(seq)+ 16)
print 'This is %d different than the dssp file, which was: ' %minus, (len(dsspAA)+16)
OHE.close()
print
print 
print "OHE found in file: ", OHE
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
out=open(sys.argv[4],'w')
np.savetxt(out, np.around(final, decimals=0), fmt='%.0f')
out.close()
print 'Sliding table found at: ', out

print 'Sliding table/ 20 should also == the length of the OHE table. Sliding table/20 =', (len(final)/20) 
