import tables
import numpy as np
import sys

feat= np.loadtxt(sys.argv[1])
OH= np.loadtxt(sys.argv[2])
#sysargv[3] is out file

index= []
for num, line in enumerate(OH):
    if num != 0 and num % 20 == 0:
        index.append(num)
OH= np.vsplit(OH, index)


h5 = tables.open_file(sys.argv[3], 'w')
one_hot = h5.create_earray(h5.root, name='one_hot', shape=(0, 20, 15), atom=tables.Int8Atom())
#pssm = h5.create_earray(h5.root, name='pssm', shape=(0, 15, 21), atom=tables.Float32Atom())
ss = h5.create_earray(h5.root, name='ss', shape=(0, 3), atom=tables.Int8Atom())

#ss.append(np.array([1, 0, 0], dtype=np.int8)[np.newaxis,:])
#one_hot.append(...)

for f,array in zip(feat,OH):
    print f
    one_hot.append(array[np.newaxis,:])
    ss.append(f[np.newaxis,:])


h5.close()






# Example of reading:

#h5 = tables.open_file('table.h5')
#print h5.root.ss[0]
#print h5.root.ss[1]
#print h5.root.ss[2]

#print
#print
#print h5.root.ss[0:2]
#print h5.root.one_hot[0:2]
