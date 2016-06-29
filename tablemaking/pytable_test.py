#1HUFA_enc is the encoded file
#simple_enc_1HUFA is the output of the sliding tables


### this is to test making a table in pytables with these two items to see if we can get what we want

import numpy as np
import sys
from tables import *

enc= open(sys.argv[1]).read().splitlines()
slide= open(sys.argv[2])  #.read().splitlines()




for n,line in enumerate(enc):
  #  print line
    pass
print n

class Practice(IsDescription):   ##
    name= StringCol(16)

h5file=open_file('practiceenc1.5h', mode= 'w', title= 'test file') ##
group= h5file.create_group('/', 'group1', 'just some group') ##
table= h5file.create_table(group, 'encoded', Practice, 'one column encoding practice') ##

value=table.row

for line in enc:
    line= line.replace(']', ',')
    line= line.replace('[', ' ')
    line= line.split(',')
    print line[1]
    if line[1]:
        value['name'] = line[1]
        value.append()
table.flush()

















'''
g= np.load(slide)
for i in g:
    x= g[i]

print x
import numpy as np
import tables

# Generate some data
#x = np.random.random((100,100,100))

# Store "x" in a chunked array...
f = tables.openFile('test.hdf', 'w')
atom = tables.Atom.from_dtype(x.dtype)
ds = f.createCArray(f.root, 'somename', atom, x.shape)
ds[:] = x
f.close()
'''
