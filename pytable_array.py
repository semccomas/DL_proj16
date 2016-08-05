import tables
import numpy as np
import sys

feat= np.loadtxt(sys.argv[1])
name= 'group_' + (sys.argv[1])[21:26]
#FOR ABOVE:::::::: [21:26] for testing data #[16:21] for when using training data
OH= np.loadtxt(sys.argv[2])
#sysargv[3] is out file

print name

index= []
for num, line in enumerate(OH):
    if num != 0 and num % 20 == 0:
        index.append(num)
OH= np.vsplit(OH, index)

#h5 = tables.open_file(sys.argv[3], 'w')
#### here i make one big table, just change out for the one above if you want lots of individual tables and for gods sake dont forget to delete the file if you run again



h5= tables.open_file('big_test_table', 'a')
group= h5.create_group('/', name, 'individual group')

one_hot = h5.create_earray(group, name='one_hot', shape=(0, 20, 15), atom=tables.Int8Atom())
#pssm = h5.create_earray(h5.root, name='pssm', shape=(0, 15, 21), atom=tables.Float32Atom())

d=feat.ndim
if d == 1:
    ss = h5.create_earray(group, name='ss', shape=(0, d), atom=tables.Int8Atom())
    feat=np.reshape(feat,(-1,1))
else:
	if len(feat[0]) == 4:
	    ss = h5.create_earray(group, name='ss', shape=(0, 4), atom=tables.Int8Atom())
	elif len(feat[0]) == 3:
	    ss = h5.create_earray(group, name='ss', shape=(0, 3), atom=tables.Int8Atom())
	elif len(feat[0]) == 2:
	    ss = h5.create_earray(group, name='ss', shape=(0, 2), atom=tables.Int8Atom())






'''
one_hot = h5.create_earray(h5.root, name='one_hot', shape=(0, 20, 15), atom=tables.Int8Atom())
#pssm = h5.create_earray(h5.root, name='pssm', shape=(0, 15, 21), atom=tables.Float32Atom())
d=feat.ndim
if d == 1:
    ss = h5.create_earray(h5.root, name='ss', shape=(0, d), atom=tables.Int8Atom())
    feat=np.reshape(feat,(-1,1))
else:
    if len(feat[0]) == 3:
        ss = h5.create_earray(h5.root, name='ss', shape=(0, 3), atom=tables.Int8Atom())
    else:
        ss = h5.create_earray(h5.root, name='ss', shape=(0, 2), atom=tables.Int8Atom())


'''





for f,array in zip(feat,OH):
    one_hot.append(array[np.newaxis,:])
    ss.append(f[np.newaxis,:])
print one_hot
print ss


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






#ss.append(np.array([1, 0, 0], dtype=np.int8)[np.newaxis,:])
#one_hot.append(...)
