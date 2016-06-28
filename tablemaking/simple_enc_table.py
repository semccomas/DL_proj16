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
#print seq

letters= 'ACDEFGHIKLMNPQRSTVWY'
#print letters


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

########### the array is, without splicing, the length of the sequence on each row, (therefore the position for each aa is the column and there are 20 rows, one per aa letter)

############################## making the sliding window ###############################################

zero= np.zeros((20,15), dtype= np.int)
zeroslow= np.concatenate((zero, array), axis= 1) 
slidearray=np.append(zeroslow, zero , axis = 1)

#we did this so that we could later have a sliding window a lot easier than trying to build one after the fact. Now the array of the sequence is surrounded by 15 zeros (and these 15 0's are 20 rows each so each row in the array starts with 20 zeros


#zeros low- sequence begins at index 16  (so char 15)
#slidearray- sequence begins at index 154 (char 154 which is included when indexing)
#i.e. the print statements below will show the very end of each array, with the very last index (or very first) is the last or first value in the sequence
#print zeroslow[:,0:16]
#print slidearray[:,154:]





for a in xrange(len(slidearray[1,:])):
  ###  print a
    ae= slidearray[:,a:a+15]
    if len(ae[1,:])== 15:
        print ae



#print len(slidearray[1,:])










### THOUGH WE DONT NEED THIS ANYMORE I am keeping it in because its a useful loop of saying, do this until that 
   # while len(ae[1,:]) != 15:
    #    ae= np.insert(ae, 0, 15, axis= 1)    #insert 15's at the 0 end of the array (left hand side) until the array is length 15
    #    if len(ae[1,:]) == 15:
      #      break
 #   print ae




