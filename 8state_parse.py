### parse 8 state, pretty much the same script as biopythondssp ++ pytable_array but for 8 states not 3. Rest is the same 

from Bio.PDB import DSSP, PDBParser 
import sys 
from difflib import * 
import numpy as np
import tables as tb

##############################################################################################################################
###################################################### PARSE DSSP ############################################################
##############################################################################################################################

fa= open(sys.argv[1]).read().splitlines()
filename=sys.argv[2]    #this will eventually become the pdb file that we run dssp on
dssp=open(filename)

p=PDBParser()
structure = p.get_structure(filename, dssp)
model= structure[0]        #there is only one structure for dssp (NMR for example has more) and the dssp parser can only take one structure
dssp= DSSP(model, filename)
a_key = list(dssp.keys())

statedic= {'H':1, 'I':2 , 'G':3, 'E':4, 'B':5, 'T':6, 'S':7, 'L':8, '-':0}
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

d=Differ()
diff= d.compare(seq, dsspAA)
comp= '\n'.join(diff)
match=0
total= [ ]

for diff in comp.split('\n'):
	if '+ X' in diff:
		match= match + 1
	elif '-' not in diff and 'X' not in diff:
		total.append(states[match])
		match= match+1
	else:
		total.append(np.nan)

padded=np.lib.pad(total, (8,8), 'constant', constant_values=(np.nan, np.nan))    #hey dipshit this is for padding the OHE features. Len of the sliding table has to == 15 always so instead of 15 on each side its 8

encoded = np.zeros((len(padded), 8))

encoded[np.where(padded == 1), 0] = 1
encoded[np.where(padded == 2), 1] = 1
encoded[np.where(padded == 3), 2] = 1
encoded[np.where(padded == 4), 3] = 1
encoded[np.where(padded == 5), 4] = 1
encoded[np.where(padded == 6), 5] = 1
encoded[np.where(padded == 7), 6] = 1
encoded[np.where(padded == 8), 7] = 1



 
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


print np.shape(final), np.shape(encoded)



##############################################################################################################################
################################################### PSSM TO PYTABLE ####################################################
##############################################################################################################################


name= 'group_' + sys.argv[1][-8:-3] 
print name 

#feature= np.loadtxt(sys.argv[2])    #this is the features for secondary structure  
h5= tb.open_file('8_state_table', 'a')    ##########!!!!!!!!!! THIS IS THE ONLY THING YOU HAVE TO CHANGE !!!!!!!!!!!!!!!! 
group= h5.create_group('/', name, 'individual group')

one_hot = h5.create_earray(group, name='one_hot', shape=(0, 20, 15), atom=tb.Float32Atom())   #would be 0, 21, 15 if you want it to be the shape of the old one
ss = h5.create_earray(group, name='ss', shape=(0, 8), atom=tb.Int8Atom())

index= []
#### splitting the sliding table into bits of 21 sized timesteps ## 
for num, line in enumerate(final):
    if num != 0 and num % 20 == 0:
        index.append(num)
final= np.vsplit(final, index)
#print np.shape(final)

for feat,line in zip(encoded, final):
    ss.append(feat[np.newaxis,:])
    one_hot.append(line[np.newaxis,:])


print ss
print one_hot

'''
h5.close()
dataset = tb.open_file('pssm_table')
for group in dataset.walk_groups():
    for array in group:
        a=array.pssm.read()
        print np.shape(a), a
'''
h5.close()