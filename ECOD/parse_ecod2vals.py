import sys
import numpy as np 
import operator 

protein_names = ((open(sys.argv[1]).read()).lower()).splitlines()
ecod = open(sys.argv[2]).read().splitlines()



#### parse each file, for prot we want the pdbID to only be 4 chars, not 5 so we take away the last to match ecod
protein = []
for line in protein_names:
	protein.append(line[0:-1])
'''
pdb= []
for line in ecod[1:]:
	line = line.split('\t')
	pdb.append(line)
pdb= np.asarray(pdb)
pdb = np.column_stack((pdb[:,4:5], pdb[:,9:10]))
'''
## above was from the old ecod files, now we will do from the new ones... Same process just different file to parse
pdb = []
for line in ecod:
		if '>' in line:
			line = line.split('|')
			dots = line[2].split('.')
			tup=[line[1][1:6], dots[0]]
			pdb.append(tup)
pdb= np.asarray(pdb)

### compare the two files, make a dict so that there are not more than one of each protein name(which there will be)
pdb_dict = {}
for line in pdb:
    if line[0][0:4] in protein:
        pdb_dict[line[0][0:4]]= line[1]

### invert dictionary above to filter by x_class name. Note that this will just sort at random or probably take the last occuring one on the list. You will have to find another way if you need to make this orderly
## as far as I can tell, the two below functions do the exact same thing in the exact same order. 
#inv_map = {v: k for k, v in pdb_dict.items()}
#x_names= dict((v, k) for k, v in pdb_dict.iteritems())

## this is where this script changes from parse_ecod... 

#first sort by class number

sorted_pdb = sorted(pdb_dict.items(), key=operator.itemgetter(1))
end_count = len(sorted_pdb)
doubles= []
for n, (prot, catg) in enumerate(sorted_pdb):
    if n +1 != end_count:
        if sorted_pdb[n][1] != sorted_pdb[n+1][1]:
            doubles.append(sorted_pdb[n])
            doubles.append(sorted_pdb[n+1])

### above will get (assuming there is more than one value per class i.e. there are 3 in class 903 or something, this will print value 1 and 3 since they are on the edges)
### below takes care of the single values...
''''
A, 901
B, 902
C, 903
D, 903
E, 903
F, 904
 
when in sorted_pdb

becomes...
... A, 901
B, 902
B, 902
C, 903
E, 903
F, 904 ...

and then the set function will remove all doubles, (i.e. B)

Essentially this is a long way of having 2 values per class, or if there is only 1 available, then just 1 value in that class
 

'''

write= list(set(doubles))
### just making it easier to parse in bash::
write = np.asarray(write, dtype = str)
for line in write:
    line[0]= line[0].upper()
np.savetxt(sys.argv[3], write, delimiter = ';' ,fmt="%s")

