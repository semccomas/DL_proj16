### this will take the .fa file and the dssp file and compare the two

import sys
from difflib import *
from sklearn.preprocessing import OneHotEncoder
import numpy as np

fa= open(sys.argv[1]).read().splitlines()
dssp=open(sys.argv[2]).read().splitlines()

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
print len(fAA), 'len faa', fAA                     #
print 'dssp', dsspAA, len(dsspAA)                  #
print x3states, len(x3states), 'states'            #
                                                   #                                
####################################################



###################################################### comparing the two files, the actual protein (fa) and the dssp file ###################

d=Differ()
diff= d.compare(fAA, dsspAA)
comp= '\n'.join(diff)

out_pre_encoded= open('not_encoded_test', 'w')                 #will become encoded numpy array later on
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



#out=open(sys.argv[3], 'w')


load=np.loadtxt(open('not_encoded_test', 'r'), dtype=str, delimiter= ',')



NumsOnly= []
for i in (load[:,2]).reshape(-1, 1):
    if '!!' not in i:
        NumsOnly.append(i)
NumsOnly= np.asarray(NumsOnly).reshape(-1,1)
print NumsOnly
        
enc= OneHotEncoder()
encoded= enc.fit_transform(NumsOnly).toarray()
encoded= np.around(encoded, decimals=0)


null= np.zeros(3)
print null

match= 0
for y,z in zip(load, encoded):
    if '!!' not in y:
        print y, encoded[match]
        match= match +1
    else:
        print y, null
    
#### I think the above zipped loop will eventually replace load with the tables and then we are donezo. 









'''
#and lastly put together the two arrays again, with the amino number on the leftmost column, and the other 3 columns are the encoded function
final=np.concatenate((amino, encoded), axis=1)


#### the fmt and the decimals are options to only print 0.00 not 0.00000000 e00 kind of thing. More manageable files


#np.savetxt(out, np.around(final, decimals=2),fmt='%.2f',delimiter='\t')

#out.close()

'''





