
import matplotlib.pyplot as plt
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten, Convolution1D, MaxPooling1D, Input, merge
#from keras.callbacks import TensorBoard
import numpy as np
import tables as tb
import sys 

#np.random.seed(3)

##############################################################################################################################
###################################################### LOADING TABLES ############################################################
##############################################################################################################################
def inputs (filename, ohe, ss):
	dataset = tb.open_file(filename)
	for group in dataset.walk_groups():
		for array in group:
			try:
				where = np.where(np.max(array.ss.read(), axis=1) == 1)
				ohe.append(array.one_hot.read()[where])		### changed from array.one_hot.read() for pssm, one more down low to change too
				ss.append(array.ss.read()[where]) 
			except AttributeError:
				pass
	dataset.close()

######################## sequence info ##########################
table= 'big_table'			#'big_table' == biopythondssp original	8state_table is 8 states
test_table= 'big_test_table'					#'big_test_table' == biopythondssp original		8_state_test_table is 8 states		
train_ohe= []
train_ss=[]
inputs(table, train_ohe, train_ss)
train_ohe= np.concatenate(train_ohe, axis = 0)			
train_ss= np.concatenate(train_ss, axis=0)
test_ohe= []
test_ss= []
inputs(test_table, test_ohe, test_ss)
test_ohe= np.concatenate(test_ohe, axis = 0)			
test_ss= np.concatenate(test_ss, axis=0)


######################## pssm info ##########################
pssm_table= 'pssm_table_jhE3'			#'big_table' == biopythondssp original			pssm_table_jhE0 pssm_table_jhE3 == pssm 
pssm_test_table= 'pssm_test_table_jhE3'					#'big_test_table' == biopythondssp original				pssm_test_table_jhE0 pssm_test_table_jhE3 == pssm
train_pssm= []
train_ss=[]
inputs(pssm_table, train_pssm, train_ss)
train_pssm= np.concatenate(train_pssm, axis = 0)			
train_ss= np.concatenate(train_ss, axis=0)
test_pssm= []
test_ss= []
inputs(pssm_test_table, test_pssm, test_ss)
test_pssm= np.concatenate(test_pssm, axis = 0)			
test_ss= np.concatenate(test_ss, axis=0)


### NOTE! Obviously the ss is the same for both pssm and seq so it does nothing to just redefine the list over again. Easier than changing the function

'''
#X_train= train_ohe  == train_pssm
#Y_train= train_ss 	 
#X_test= test_ohe	== test_pssm
#Y_test= test_ss

print np.shape(train_pssm), 'pssm'
print np.shape(train_ss), 'sstr'
print np.shape(train_ohe), 'ohe'
print np.shape(test_ss), 'ss'
print np.shape(test_pssm), 'tespss'
print np.shape(test_ohe), 'testohs'
'''
print 
##############################################################################################################################
###################################################### THE MODEL ############################################################
##############################################################################################################################
seq_input = Input(shape=(20,15))
cs = Convolution1D(10, 3, activation= 'relu')(seq_input)
cs = Convolution1D(10, 3, activation= 'relu')(cs)
#cs = Convolution1D(10, 3, activation= 'relu')(cs)
cs = Convolution1D(10, 3, activation= 'relu')(cs)
#cs = Dropout(0.5)(cs)
sx = Flatten()(cs)

pssm_input = Input(shape=(21,15)) 
cp = Convolution1D(10, 3, activation= 'relu')(pssm_input) 
cp = Convolution1D(10, 3, activation= 'relu')(cp)
#cp = Convolution1D(10, 3, activation= 'relu')(cp)
cp = Convolution1D(10, 3, activation= 'relu')(cp)
px = Flatten()(cp)

x = merge([sx, px], mode='concat')
x = Dense(64, activation='relu')(x)
x = Dense(64, activation='relu')(x)
x = Dropout(3)(x)
predictions = Dense(3, activation='softmax')(x)

model = Model(input=[seq_input, pssm_input], output=predictions)		
model.compile(optimizer='adagrad',
              loss='categorical_crossentropy',
              metrics=['accuracy'])
model.summary()
nb_epoch= 20
batchsize= 200

history = model.fit([train_ohe, train_pssm], train_ss, nb_epoch=nb_epoch, batch_size=batchsize, validation_data=([test_ohe, test_pssm], test_ss))


from keras.utils.visualize_util import plot
plot(model, to_file='model.png')




loss=history.history['loss']    
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

#scores = model.evaluate(test_pssm, test_ss, batch_size=batchsize)
#print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100)), 'accuracy'    #prints accuracy, one value, in percent 

##############################################################################################################################
################################################### EVALUATION / PLOTTING ####################################################
##############################################################################################################################

#### plotting loss ####
epoch= []
for r in xrange(nb_epoch):
    epoch.append(r)

plt.figure(3)
plt.plot(epoch, loss)
plt.ylabel('loss')
plt.xlabel('epoch')
plt.ylim(ymin=0)

plt.show()


pssmdata = tb.open_file(pssm_test_table)
seqdata = tb.open_file(test_table)
for pgroup , sgroup in zip(pssmdata.walk_groups() , seqdata.walk_groups()):
	for parray, sarray in zip (pgroup, sgroup):
		pssm = parray.one_hot.read()
		seq = sarray.one_hot.read()
		#single_protein= array.one_hot.read()		## this I also changed from array.one_hot.read() to pssm. Should only be two instances in whole file
		#print np.shape(single_protein)
		predict=model.predict([seq, pssm], batch_size= batchsize, verbose= 1)
		raw_input()

		maximum= np.amax(predict, axis= 1)				# max probability in each row
		positions= np.argmax(predict, axis= 1)			# which index that (^) probability is. Is same as model.predict_classes
		

		# compares two lists, list 2 is an index for the 3 values in list one, appends the maximum value to the appropriate group. Rest is 0. 
		# pos is to later use in graph to show if predicted to be present. Really the only diff is threshold
		def find (lst, pos, number):
			for index in xrange(len(positions)):
				if positions[index] == number:
					pos.append(0.5)					#0.5 is here just because its in the middle of the graph so you can see it better. Can be whatever
		 			if maximum[index] >= 0.6:		#### !!!!!!!!!!!!! This is the one you can change if you want it to be a different threshold
						lst.append(maximum[index])		
					else:
						lst.append(0)
				else: 
					lst.append(0)
					pos.append(np.nan)         #nans dont get plotted so this removes the part of the plot with no info


		helix= []
		h_pos= []
		sheet= []
		s_pos= []
		coil= []
		c_pos=[]		
		find(helix, h_pos, 0)
		find(sheet, s_pos, 1)
		find(coil, c_pos, 2)

		
		print len(sheet), len(s_pos)
		print len(helix), len(h_pos)
		print len(coil), len(c_pos)
	
		f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
		ax1.plot(helix, 'r')
		ax1.plot(h_pos, 'ko')
		ax1.set_title('Helix')
		ax2.plot(sheet, 'g')
		ax2.plot(s_pos, 'ko')
		ax2.set_title('Sheet')
		ax3.plot(coil, 'b')
		ax3.plot(c_pos, 'ko')
		ax3.set_title('Coil')

		plt.ylim(ymin= -1, ymax= 2)
		plt.xlabel('Position')
		plt.ylabel('Probability' )
		plt.legend(loc= 'upper left')
		
		plt.show()

		#plt.savefig()

		################ index 0 is helix, index 1 is strand and index 2 is coil 
		####'H':0, 'I':0 , 'G':0, 'E':1, 'B':1, 'T':2, 'S':2, 'L':2
		name= parray._v_name
		print name

pssmdata.close()
seqdata.close()
		
