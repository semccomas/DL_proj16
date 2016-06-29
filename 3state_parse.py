### this will take the .fa file and the dssp file and compare the two

import sys
from difflib import *

fa= open(sys.argv[1]).read().splitlines()
dssp=open(sys.argv[2]).read().splitlines()

aminodic= {'A':1, 'R':2, 'N':3, 'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20} 
statedic={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2}
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






##### just for comparing in writing the program #####
fAA= fa[1]                                         #
print len(fAA), 'len faa', fAA                     #
print 'dssp', dsspAA, len(dsspAA)                  #
print states, len(states), 'states'                #
                                                   #                                #
####################################################



###################################################### comparing the two files, the actual protein (fa) and the dssp file ###################

d=Differ()
diff= d.compare(fAA, dsspAA)
comp= '\n'.join(diff)


match= 0                                       #need this here as well as enumerate because we only want to count matches
for n, line in enumerate(comp.splitlines()):
    if '-' not in line:
        print line[2], dsspAA[match], states[match]
        match= match +1 
    else:
        print fAA[n], '!!', '!!'
print 'line[2], dsspAA[match], states[match]'

#here !! symbolizes blank or missing, for both the amino acid and the states

######################################## one hot encoding the states #####################################

