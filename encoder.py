
#### this takes the parsed but not encoded (dssp_values_not_encoded) as an input 1 and the file to write with as input 2 (dssp_values) and transforms this into a one hot encoded array

import sys
from sklearn.preprocessing import OneHotEncoder
import numpy as np

#sys.argv[1]= input to array
out=open(sys.argv[2], 'w')

### load the text as an array, make the arrays with the same shape (will look like: 
#[[0],
#[1],
#[1], .... 

array=np.loadtxt(sys.argv[1], delimiter= ',')
states= (array[:,1]).reshape(-1, 1)   #'reshape your array if it has one feature to -1, 1 or 1,-1 if it has one sample.'......
# Would have thought this is only one sample but this seems to be what I want for the time being
amino=(array[:,0]).reshape(-1,1)

#transform the 3 states to one hot encoded
enc= OneHotEncoder()
encoded= enc.fit_transform(states).toarray()
encoded= np.around(encoded, decimals=3)
print encoded

#and lastly put together the two arrays again, with the amino number on the leftmost column, and the other 3 columns are the encoded function
final=np.concatenate((amino, encoded), axis=1)
#### the fmt and the decimals are options to only print 0.00 not 0.00000000 e00 kind of thing. More manageable files
np.savetxt(out, np.around(final, decimals=2),fmt='%.2f',delimiter='\t')

out.close()
