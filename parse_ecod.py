import sys
import numpy as np 

protein_names = ((open(sys.argv[1]).read()).lower()).splitlines()
ecod = open(sys.argv[2]).read().splitlines()



#### parse each file, for prot we want the pdbID to only be 4 chars, not 5 so we take away the last to match ecod
protein = []
for line in protein_names:
	protein.append(line[0:-1])

pdb= []
for line in ecod[5:]:
	line = line.split('\t')
	pdb.append(line)
pdb= np.asarray(pdb)
pdb = np.column_stack((pdb[:,4:5], pdb[:,9:10]))



### compare the two files, make a dict so that there are not more than one of each protein name(which there will be)
pdb_dict = {}
for line in pdb:
    if line[0] in protein:
        pdb_dict[line[0]]= line[1]

### invert dictionary above to filter by x_class name. Note that this will just sort at random or probably take the last occuring one on the list. You will have to find another way if you need to make this orderly
## as far as I can tell, the two below functions do the exact same thing in the exact same order. 
#inv_map = {v: k for k, v in pdb_dict.items()}
x_names= dict((v, k) for k, v in pdb_dict.iteritems())
write = []
for key, val in x_names.iteritems():
	tup = [key.upper, val.upper()]
	write.append(tup)
### just making it easier to parse in bash::
write = np.asarray(write, dtype = str)
np.savetxt(sys.argv[3], write, delimiter = ';' ,fmt="%s")

#out= open(sys.argv[3], 'w')
#out.write(x_names)
#out.close()