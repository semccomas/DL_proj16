######################## THIS IS THE BIG GUY, DOES BOTH ONE HOT ENCODING (FROM 3STATEPARSE) AND THE SLIDING AMINO ACID PROFILE (FROM SLIDING_TABLE)

### this will take the .fa file and the dssp file and compare the two. There are two outputs, the one hot encoded profile for the states and the table

import sys
from difflib import *
from sklearn.preprocessing import OneHotEncoder
import numpy as np

fa= open(sys.argv[1]).read().splitlines()
dssp=open(sys.argv[2]).read().splitlines()

## LEAVE THESE TWO COMMENTED OUT, JUST FOR REFERENCE TO SEE WHAT WE USE: encoded_outfile= open(sys.argv[3], 'w')
## out= open(sys.argv[4], 'w')


aminodic= {'A':1, 'R':2, 'N':3, 'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20} 
statedic={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2, '!!': '!!'}
#0 is helix, 1 is strand, 2 is coil

#################################################### parse the DSSP output and find the states and the amino acid ###############################################            
dsspAA= []
states= []
start=0
for line in dssp:
    value= line.split()
    if 'N-H-->O' in value:                       #indicates where to start parsing
        start= 1 
        continue                              #skip the header line
    if not start: 
        continue                              #skip all lines before the header line also



    if len(value[1])==4 and value[2] in aminodic:
        dsspAA.append(value[2])
        if value[3] in statedic:
            states.append(value[3])
        elif value[3] not in statedic:
            states.append('!!') #### same as below, no state, how to encode??
    elif len(value[1])!=4 and value[3] in aminodic:
        dsspAA.append(value[3])
        if value[4] in statedic:
            states.append(value[4])
        elif value[4] not in statedic:
            states.append('!!') ### this is if there is no state but there is an aa. what to do to be encoded 000??



###### transform the states to only 3 instead of 8. 
#In the future you can probably just take this away and do one hot encoding directly
x3states= []
for s in states:
    x3states.append(statedic[s])



##### just for comparing in writing the program #####
fAA= fa[1]                                         #
#print len(fAA), 'len faa', fAA                     #
#print 'dssp', dsspAA, len(dsspAA)                  #
#print x3states, len(x3states), 'states'            #
                                                   #                                
####################################################



###################################################### comparing the two files, the actual protein (fa) and the dssp file ###################

d=Differ()
diff= d.compare(fAA, dsspAA)
comp= '\n'.join(diff)

out_pre_encoded= open('not_encoded_forloadingonly', 'w')                 #will become encoded numpy array later on
match= 0                                                       #need both a counter as well as enumerate to be accessing two different lists with diff. values
for n, line in enumerate(comp.splitlines()):
    if '-' not in line:
        out_pre_encoded.write(line[2] + ',' + dsspAA[match] + ',' +  str(x3states[match]) + '\n')
       # print line[2], dsspAA[match], x3states[match]
        match= match +1 
    else:
        out_pre_encoded.write(fAA[n] + ',' + '!!' + ',' + '!!' + '\n')
       # print fAA[n], '!!', '!!'

#here !! symbolizes blank or missing, for both the amino acid and the states

out_pre_encoded.close()


######################################## one hot encoding the states ###################################################################################

#### this takes the parsed but not encoded (dssp_values_not_encoded) as an input 1 and the file to write with as input 3 (dssp_values) and transforms this into a one hot encoded array........... this is also modified from encoder.py

### load the text as an array, make the arrays with the same shape (will look like: 
#[[0],
#[1],
#[1], ....reshape your array if it has one feature to -1, 1 or 1,-1 if it has one sample.. We want -1, 1
#######################################################################################################################################################



load=np.loadtxt(open('not_encoded_forloadingonly', 'r'), dtype=str, delimiter= ',')

NumsOnly= []
for i in (load[:,2]).reshape(-1, 1):
    if '!!' not in i:
        NumsOnly.append(i)
NumsOnly= np.asarray(NumsOnly).reshape(-1,1)
#print NumsOnly
        
enc= OneHotEncoder()
encoded= enc.fit_transform(NumsOnly).toarray()
encoded= np.around(encoded, decimals=0)


null= np.zeros(3)

encoded_outfile= open(sys.argv[3], 'w')
match= 0
for y,z in zip(load, encoded):
    if '!!' not in y:
        #print y, encoded[match]
        encoded_outfile.write(str(y) + ' ' + str(encoded[match]) + '\n')
        match= match +1
    else:
        #print y, null
        encoded_outfile.write(str(y) + ' ' + str(null) + '\n')

encoded_outfile.close()
 

print 'One hot encoding found in file: ' , encoded_outfile   
 
##############################################################################################################################
###################################################### SLIDING TABLE #########################################################
##############################################################################################################################

############### RUN THIS ON sysargv1


######################### importing the sequence and extracting information #############################
### seq is the actual sequence
### letters are all 20 amino acid letters

#f= open(sys.argv[1]).read().splitlines()
out= open(sys.argv[4], 'w')
for line in fa:
    if '>' not in line:
        seq= line

letters= 'ACDEFGHIKLMNPQRSTVWY'


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




#out= open(sys.argv[2], 'w')

vstack= []

for a in xrange(len(slidearray[1,:])):
  ###  print a
    ae= slidearray[:,a:a+15]
    if len(ae[1,:])== 15:
        vstack.append(ae)
        #print ae
       # np.savez(out, ae) 



vSTACK= np.vstack(vstack)
np.savetxt(out, np.around(vSTACK, decimals=0), fmt='%.0f')
#print vSTACK

out.close()

##print len(slidearray[1,:])




print 'Sliding table found at: ',  out

