from Bio.PDB import DSSP, PDBParser 
import sys 
from difflib import * 
import numpy as np
import tables as tb
import math
import collections

##############################################################################################################################
###################################################### PARSE DSSP ############################################################
##############################################################################################################################

fa= open(sys.argv[1]).read().splitlines()
filename=sys.argv[2]    #this will eventually become the pdb file that we run dssp on
#### alignment is sys.argv[3]

dssp=open(filename)

p=PDBParser()
structure = p.get_structure(filename, dssp)
model= structure[0]        #there is only one structure for dssp (NMR for example has more) and the dssp parser can only take one structure
dssp= DSSP(model, filename)
a_key = list(dssp.keys())

statedic3= {'H':1, 'I':2 , 'G':3 , 'E':4 , 'B':5, 'T':6, 'S':7, '-':8} #, '-':0}
statedic8={'H':60, 'I':60, 'G':60, 'E':80, 'B':80, 'T':1, 'S':1, '-':1}
dsspAA=[]
states3= [] 
states8= []
rsa= []
for line in a_key:
    dsspAA.append(dssp[line][1])
    states3.append(statedic3[dssp[line][2]])
    states8.append(statedic8[dssp[line][2]])
    rsa.append(dssp[line][3])


rsa= np.asarray(rsa)
rsa[np.where(rsa == 'NA')] = np.nan
rsa=np.asarray(rsa, dtype= float)


##############################################################################################################################
###################################################### PARSE PSSM #######################################################
##############################################################################################################################


def parse_PSSM(alignment):
    one2number = 'ARNDCEQGHILKMFPSTWYV-'
    bi = [0.0825, 0.0553, 0.0406, 0.0545, 0.0137, 0.0393, 0.0675, 0.0707, 0.0227, 0.0595, 0.0966, 0.0584, 0.0242, 0.0386, 0.0470, 0.0657, 0.0534, 0.0108, 0.0292, 0.0687]
    b = dict()
    for index, letter in enumerate(one2number[:-1]):
        b[letter] = bi[index]

    counts = collections.defaultdict(lambda: collections.defaultdict(int))
    seqcount = 0.
    gapcount = 0
    
    for line in open(alignment):
        line = line.strip()
        seqcount += 1
        for position, letter in enumerate(line):
            counts[position][letter] += 1
        gapcount += line.count('-')
        line_length = len(line)

    b['-'] = gapcount/(seqcount * line_length)

    pssm = np.zeros((len(counts.keys()), len(one2number)), dtype=np.float32)
    pssm2 = np.zeros((len(counts.keys()), len(one2number)), dtype=np.float32)
    for position in counts:
        q0 = []
        q1 = []
        for letter in one2number:
            p = counts[position][letter] / (b[letter] * seqcount)
            if p > 0:
                q0.append(math.log(p))
                q1.append(p * math.log(p))
            else:
                q0.append(math.log(0.1 / (b[letter] * seqcount)))
                q1.append(0)

        pssm[position, :] = q0
        pssm2[position, :] = q1

    return pssm #,  pssm2

######## This array was originally shaped as (seqlen, 21). Reshapign below leads it to be: (21, seqlen), which was the same shape
#as the array in biopythondssp.py before making it into a sliding table. Below I have taken parts of the code from the sliding table in order to make
## it the same as the biopythondssp output in hopes of giving myself less of a headache later. It should afterwards be the shape: (seqlen, 21, 15) where
# 21 is the timestep, i.e. the number of features (aa alphabet plus '-') and 15 is the len of the window.
#Only real difference than the old file from biopythondssp is that instead of all 0's except a 1, these are all the PSSM values

## make the array like the old one, and pad it with zeros on both sides, seq length is therefore 30 longer, make it sliding. Should be exact same shape as before
pssm= np.fliplr(np.rot90(parse_PSSM(sys.argv[3]), 3))
seq = fa[1]

##### encoding the sequence data to be a table
letters= 'ACDEFGHIKLMNPQRSTVWY'   # to create numbers for each amino acid
def encode (seq, l):
    row = []
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

##### for both sequence info and pssm (WholeSeq and pssm)
def sliding_table(array):
	padded_array= []
	for line in array:
	    pad=np.lib.pad(line, (15,15), 'constant', constant_values=(0,0))
	    padded_array.append(pad)
	padded_array= np.asarray(padded_array)

	window=[]
	for index in xrange(np.shape(padded_array)[1]):       # should be len of pssm , or WholeSeq + 30!!!!!! 
	        slide= padded_array[0:,index:index+15]               #0: == use all rows ',' <column no> ':' <column no>
	        if len(slide[0]) == 15:
	            pass
	            for line in slide:                  # because you cant append arrays to each other, change to a list and then back to array. Doesn't seem to hurt the values at all
	                window.append(line)
	window=np.asarray(window)
	return window

WholeSeq = sliding_table(WholeSeq)
pssm = sliding_table(pssm)

print np.shape(WholeSeq)
print np.shape(pssm)


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

OHE=open(sys.argv[3], 'w')
np.savetxt(OHE, np.around(encoded, decimals=0), fmt='%.0f')
OHE.close()




##############################################################################################################################
################################################### PSSM TO PYTABLE ####################################################
##############################################################################################################################


name= 'group_' + sys.argv[1][-8:-3] 
print name 

#feature= np.loadtxt(sys.argv[2])    #this is the features for secondary structure  
h5= tb.open_file('8state_table', 'a')    ##########!!!!!!!!!! THIS IS THE ONLY THING YOU HAVE TO CHANGE !!!!!!!!!!!!!!!! 
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