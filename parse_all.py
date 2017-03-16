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

statedic3={'H':1, 'I':1, 'G':1, 'E':2, 'B':2, 'T':3, 'S':3, '-':3}
statedic8= {'H':1, 'I':2 , 'G':3 , 'E':4 , 'B':5, 'T':6, 'S':7, '-':8} #, '-':0}
dsspAA=[]
states3= [] 
states8= []
rsa= []
for line in a_key:
 #   print dssp[line]
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

#############################################################################################################################
###################################################### SLIDING TABLES #######################################################
##############################################################################################################################

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

print np.shape(WholeSeq), 'shapeWS'
print np.shape(pssm), 'shapePSSM'

##############################################################################################################################
######################################################ONE HOT ENCODING #######################################################
##############################################################################################################################

### run this function for RSA, 3 state, and 8 state!
# the constant can just stay nan. For RSA, we will keep this nan and use it later, for the one hot encoding there will be no nan because it will be replaced with zeros!!!! 
def compare_files(feature): 		
	d=Differ()
	diff= d.compare(seq, dsspAA)
	comp= '\n'.join(diff)
	match=0
	total= [ ]

	for diff in comp.split('\n'):
		if '+ X' in diff:
			match= match + 1
		elif '-' not in diff and 'X' not in diff:
			total.append(feature[match])
			match= match+1
		else:
			total.append(np.nan)

	total=np.lib.pad(total, (8,8), 'constant', constant_values=(np.nan, np.nan))    #this is for padding the OHE features. Len of the sliding table has to == 15 always so instead of 15 on each side its 8
	return total
	#return comp

states3 = compare_files(states3)
states8 = compare_files(states8)
rsa = compare_files(rsa)
#print rsa

#### these 3 should all be the same shape. That is, pssm/ WholeSeq divided by 21 or 20, respectively. Also = sequence length + 16
### this below only pertains to 3 and 8 states. The output of compare_files(total) is the input to this

def OHE (array, nb_feats):
	encoded = np.zeros((len(array), nb_feats))

	encoded[np.where(array == 1), 0] = 1
	encoded[np.where(array == 2), 1] = 1
	encoded[np.where(array == 3), 2] = 1
	if nb_feats == 8:
		encoded[np.where(array == 4), 3] = 1
		encoded[np.where(array == 5), 4] = 1
		encoded[np.where(array == 6), 5] = 1
		encoded[np.where(array == 7), 6] = 1
		encoded[np.where(array == 8), 7] = 1

	return encoded 

states3 = OHE(states3, 3)
states8 = OHE(states8, 8)


##############################################################################################################################
################################################### TO PYTABLE #############################################################
##############################################################################################################################


name= 'group_' + sys.argv[1][-8:-3] 
print name 

#feature= np.loadtxt(sys.argv[2])    #this is the features for secondary structure  
h5= tb.open_file(sys.argv[4], 'a')
group= h5.create_group('/', name, 'individual group')

seq_tab = h5.create_earray(group, name='seq_tab', shape=(0, 20, 15), atom=tb.Float32Atom())   #would be 0, 21, 15 if you want it to be the shape of the old one
pssm_tab = h5.create_earray(group, name='pssm_tab', shape=(0, 21, 15), atom=tb.Float32Atom())   #would be 0, 21, 15 if you want it to be the shape of the old one
ss3_feat = h5.create_earray(group, name='ss3_feat', shape=(0, 3), atom=tb.Int8Atom())
ss8_feat = h5.create_earray(group, name='ss8_feat', shape=(0, 8), atom=tb.Int8Atom())
rsa_feat = h5.create_earray(group, name='rsa_feat', shape=(0, 1), atom=tb.Float32Atom())
rsa=np.reshape(rsa,(-1,1))


############## might need to add filter!!!!! ################# 
### also there were 3 difficult files, something went wrong in comp. In the future we can maybe take these away with an if loop just to check that they match




#### splitting the sliding table into bits of 21 (pssm) or 20 (seq) sized timesteps ## 
def timesteps (array, size):
	index= []
	for num, line in enumerate(array):
	    if num != 0 and num % size == 0:
	        index.append(num)
	array= np.vsplit(array, index)     #split the array each time the window ends, this creates timesteps. index is the index where the window ends
	return array

WholeSeq = timesteps(WholeSeq, 20)
pssm= timesteps(pssm, 21)




##### putting all values into the pytable #############
def to_table (array, table):
	for line in array:
		table.append(line[np.newaxis,:])
	return table 

WholeSeq = to_table(WholeSeq, seq_tab)
pssm = to_table(pssm, pssm_tab)
states3 = to_table (states3, ss3_feat)
states8 = to_table (states8, ss8_feat)
rsa = to_table (rsa, rsa_feat)

