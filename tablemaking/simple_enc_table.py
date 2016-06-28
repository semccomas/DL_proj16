
############ imports 
import numpy as np
import sys


######################### importing the sequence and extracting information #############################
### seq is the actual sequence
### letters are all 20 amino acid letters

f= open(sys.argv[1]).read().splitlines()

for line in f:
    if '>' not in line:
        seq= line
print seq

letters= 'ACDEFGHIKLMNPQRSTVWY'
print letters


###################### corresponding position to amino acid and forming a table ###############################

def encode (seq, l):
    lists= []
    for c in seq:
        if c == l:
            lists.append(1)
        else:
            lists.append(0)
    return lists

array= []
for l in letters:
    array.append(encode(seq, l))
array= np.asarray(array)


############################## making the sliding window ###############################################

for a in xrange(len(seq)):
    print a
    ae= array[:,a:a+15]
    while len(ae[1,:]) != 15:
        ae= np.insert(ae, 0, 0, axis= 1)
        if len(ae[1,:]) == 15:
            break
    print ae

zero= np.zeros((20,15), dtype=np.int)
print zero
LowEnd= np.concatenate((zero, array), axis=1)
print LowEnd
np.savetxt('testfile', LowEnd)


#elif a + 15 > len(seq):
 #       ae= array[:,a:a+15]
  #      ae= np.insert(ae, 0, 3, axis= 1)



#    if len(ae[1,:]) != 15:


print len(seq)
