import sys

f= open(sys.argv[1]).read().splitlines()
start=0

amino= {'A':1, 'R':2, 'N':3, 'D':4, 'C':5, 'Q':6, 'E':7, 'G':8, 'H':9, 'I':10, 'L':11, 'K':12, 'M':13, 'F':14, 'P':15, 'S':16, 'T':17, 'W':18, 'Y':19, 'V':20, '0':0, 'X':0} 
  #previous indices for the letter to number amino acids. Here in DSSP there are sometimes blanks or '!' as the missing AA, this reads them as a 0 for some reason so I also give them the value 0 (last value in amino dict). Same with X, is nonvalid

state={'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2}
#0 is helix, 1 is strand, 2 is coil

for line in f:
    sl= line.split()
    if 'N-H-->O' in sl:                       #indicates where to start parsing
        start= 1 
        continue                              #skip the header line
    if not start: 
        continue                              #skip all lines before the header line also
   
### printing the amino accession number along with the designated state number. There are some situations where words get clumped together (if loop #2, if sl[3], so split doesnt work and we have to index earlier. 
#### basically there are 4 scenarios, where both are not right, one is not right, the other is not right, or both are right (right= at the right index and with the expected value. This accounts for all of them. I dont think order matters at all since its just a number

    if sl[3] not in amino and sl[3] not in state:
        print amino[sl[2]], ',' , 2

    elif sl[3] not in amino and sl[3] in state:    #you have to backtrack since the index is off, pretend for this one that the 3rd index is the structure state. Yes i tested this and it was right. keep sl[3] for the first two if's 
        print amino[sl[2]], ',' , state[sl[3]]

    elif sl[3] in amino and sl[4] not in state:
        print amino[sl[3]], ',' , 2

    else:
        print amino[sl[3]], ',' , state[sl[4]]

#### THIS EVEN PRINTS IT IN THE SAME ORDER AS FROM DSSP WOOP
