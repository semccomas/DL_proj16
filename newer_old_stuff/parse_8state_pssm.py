# this parser comes from David for parsing the .fa.jE??.aln files 
import math
import collections
import tables as tb 
import numpy as np
import sys

#### This is pretty much the same as parse_pssm. Just added for 8 states

##############################################################################################################################
################################################### PARSE AND MAKE PSSM ####################################################
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
seq= np.fliplr(np.rot90(parse_PSSM(sys.argv[1]), 3))

WholeSeq= []
for line in seq:
    padded=np.lib.pad(line, (15,15), 'constant', constant_values=(0,0))
    WholeSeq.append(padded)
WholeSeq= np.asarray(WholeSeq)

window=[]
for index in xrange(len(padded)):
        slide= WholeSeq[0:,index:index+15]               #0: == use all rows ',' <column no> ':' <column no>
        if len(slide[0]) == 15:
            pass
            for line in slide:                  # because you cant append arrays to each other, change to a list and then back to array. Doesn't seem to hurt the values at all
                window.append(line)

window=np.asarray(window)
#np.savetxt('pssm_parsed_%s' % sys.argv[1][16:25], array, fmt= '%d')



##############################################################################################################################
################################################### PSSM TO PYTABLE ####################################################
##############################################################################################################################


name= 'group_' + sys.argv[1][-17:-12] 

feature= np.loadtxt(sys.argv[2])    #this is the features for secondary structure  
h5= tb.open_file('pssm_test_8state_jhE0', 'a')    ##########!!!!!!!!!! THIS IS THE ONLY THING YOU HAVE TO CHANGE !!!!!!!!!!!!!!!! 
########## 4 tables here.... pssm_8state_jhE0, pssm_test_8state_jhE0, pssm_8state_jhE3, pssm_test_8state_jhE3
group= h5.create_group('/', name, 'individual group')


pssm = h5.create_earray(group, name='one_hot', shape=(0, 21, 15), atom=tb.Float32Atom())   #would be 0, 21, 15 if you want it to be the shape of the old one
ss = h5.create_earray(group, name='ss', shape=(0, 8), atom=tb.Int8Atom())

index= []
#### splitting the sliding table into bits of 21 sized timesteps ## 
for num, line in enumerate(window):
    if num != 0 and num % 21 == 0:
        index.append(num)
window= np.vsplit(window, index)
#print np.shape(window)

for feat,line in zip(feature, window):
    ss.append(feat[np.newaxis,:])
    pssm.append(line[np.newaxis,:])


print ss
print pssm

'''
h5.close()
dataset = tb.open_file('pssm_table')
for group in dataset.walk_groups():
    for array in group:
        a=array.pssm.read()
        print np.shape(a), a
'''
h5.close()