
##############################################################################################################################
###################################################### PARSE DSSP ############################################################
##############################################################################################################################
### taken from 8state_parse.py

from Bio.PDB import DSSP, PDBParser 
import sys 
import numpy as np
from difflib import *
import tables as tb

fa= open(sys.argv[1]).read().splitlines()
filename=sys.argv[2]    #this will eventually become the pdb file that we run dssp on
dssp=open(filename)

#max_val= {'A':106, 'C':135, 'D':163, 'E':194, 'F':197, 'G':84, 'H':184, 'I':169, 'K':205, 'L':164, 'M':188, 'N':157, 'P':136,  'R':248, \
#'S':130, 'T':142, 'V':142, 'W':227, 'Y':222, 'Z':196}

#'B':160,    'Q':194,      'X':222 

p=PDBParser()
structure = p.get_structure(filename, dssp)
model= structure[0]        #there is only one structure for dssp (NMR for example has more) and the dssp parser can only take one structure
dssp= DSSP(model, filename)
a_key = list(dssp.keys())


rsa= []
dsspAA= []
for line in a_key:
	#print dssp[line]
	rsa.append(dssp[line][3])
	dsspAA.append(dssp[line][1])

rsa= np.asarray(rsa)
rsa[np.where(rsa == 'NA')] = np.nan
rsa=np.asarray(rsa, dtype= float)


##############################################################################################################################
######################################################MAKING PADDED ARRAY #######################################################
##############################################################################################################################

seq= fa[1]

d=Differ()
diff= d.compare(seq, dsspAA)
comp= '\n'.join(diff)
match=0
total= [ ]

for diff in comp.split('\n'):
	if '+ X' in diff:
		match= match + 1
	elif '-' not in diff and 'X' not in diff:
		total.append(rsa[match])
		match= match+1
	else:
		total.append(np.nan)


padded=np.lib.pad(total, (8,8), 'constant', constant_values=(np.nan, np.nan))    #hey dipshit this is for padding the OHE features. Len of the sliding table has to == 15 always so instead of 15 on each side its 8

np.savetxt(sys.argv[3],padded)

 
##############################################################################################################################
###################################################### SLIDING TABLE #########################################################
##############################################################################################################################

#seq= fa[1] (just a reminder what seq is- appears in OHE)
zeros= 15        #change this accordingly to how many zeros you want on each side. Be sure to think of how many lines of zeros you will have in the OHE part
seq= ('o' * zeros) + seq + ('o' * zeros)


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
#			print slide, np.shape(slide)
			for line in slide: 					# because you cant append arrays to each other, change to a list and then back to array. Doesn't seem to hurt the values at all
				final.append(line)
final=np.asarray(final)

print np.shape(final), np.shape(padded)


##############################################################################################################################
################################################### PSSM TO PYTABLE ####################################################
##############################################################################################################################


name= 'group_' + sys.argv[1][-8:-3] 
print name 

#feature= np.loadtxt(sys.argv[2])    #this is the features for secondary structure  
h5= tb.open_file(sys.argv[4], 'a')    ##########!!!!!!!!!! THIS IS THE ONLY THING YOU HAVE TO CHANGE !!!!!!!!!!!!!!!! 
group= h5.create_group('/', name, 'individual group')

one_hot = h5.create_earray(group, name='one_hot', shape=(0, 20, 15), atom=tb.Float32Atom())   #would be 0, 21, 15 if you want it to be the shape of the old one
ss = h5.create_earray(group, name='ss', shape=(0, 1), atom=tb.Int8Atom())
padded=np.reshape(padded,(-1,1))

#### splitting the sliding table into bits of 20 sized timesteps ## 
index= []
for num, line in enumerate(final):
    if num != 0 and num % 20 == 0:
        index.append(num)
final= np.vsplit(final, index)


for feat,line in zip(padded, final):
    ss.append(feat[np.newaxis,:])
    one_hot.append(line[np.newaxis,:])


print ss
print one_hot


