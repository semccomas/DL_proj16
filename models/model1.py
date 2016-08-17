import sys 
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten, Convolution1D, MaxPooling1D, Input, merge
import numpy as np
import tables as tb
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt


##############################################################################################################################
###################################################### LOADING TABLES ############################################################
##############################################################################################################################
#def inputs (filename, seq, pssm, ss3, ss8, rsa):
#seq =[]
#pssm =[]
#ss3=[]
#ss8=[]
#rsa=[]
def inputs (filename, seq, pssm, ss3, ss8, rsa):
	dataset = tb.open_file(filename)
	for group in dataset.walk_groups():
			for array in group:
				try:
					where = np.where(np.max(array.ss3_feat.read(), axis=1) == 1)
					seq.append(array.seq_tab.read()[where])		
					pssm.append(array.pssm_tab.read()[where]) 
					ss3.append(array.ss3_feat.read()[where]) 
					ss8.append(array.ss8_feat.read()[where]) 
					rsa.append(array.rsa_feat.read()[where]) 
				except AttributeError:
					pass
	dataset.close()
#return where

######################## sequence info ##########################
table= 'table_all_jhE0'			
train_seq = []
train_pssm = []
train_ss3 = []
train_ss8 = []
train_rsa = []
inputs(table, train_seq, train_pssm, train_ss3, train_ss8, train_rsa)
train_seq= np.concatenate(train_seq, axis = 0)
train_pssm= np.concatenate(train_pssm, axis = 0)			
train_ss3= np.concatenate(train_ss3, axis = 0)			
train_ss8= np.concatenate(train_ss8, axis = 0)			
train_rsa= np.concatenate(train_rsa, axis = 0)			
	
test_table= 'test_table_all_jhE0'						
test_seq = []
test_pssm = []
test_ss3 = []
test_ss8 = []
test_rsa = []
inputs(test_table, test_seq, test_pssm, test_ss3, test_ss8, test_rsa)
test_seq= np.concatenate(test_seq, axis = 0)
test_pssm= np.concatenate(test_pssm, axis = 0)			
test_ss3= np.concatenate(test_ss3, axis = 0)			
test_ss8= np.concatenate(test_ss8, axis = 0)			
test_rsa= np.concatenate(test_rsa, axis = 0)	






##############################################################################################################################
###################################################### THE MODEL ############################################################
##############################################################################################################################
seq_input = Input(shape=(20,15))
cs1 = Convolution1D(10, 3, activation= 'relu')(seq_input)
cs2 = Convolution1D(10, 3, activation= 'relu')(seq_input)
cs3 = Convolution1D(10, 3, activation= 'relu')(seq_input)
cs4 = Convolution1D(10, 3, activation= 'relu')(seq_input)
#cs = Dropout(0.5)(cs)
cm1 = merge([cs1,cs2] , mode= 'concat')
cm2= merge([cs3,cs4], mode= 'concat')
cm= merge([cm1, cm2], mode = 'concat')
cm = Dropout(0.1)(cm)
sx = Flatten()(cm)

pssm_input = Input(shape=(21,15)) 
cp1 = Convolution1D(10, 3, activation= 'relu')(pssm_input) 
cp2 = Convolution1D(10, 3, activation= 'relu')(pssm_input)
cp3 = Convolution1D(10, 3, activation= 'relu')(pssm_input)
cp4 = Convolution1D(10, 3, activation= 'relu')(pssm_input)
pm1= merge([cp1, cp2], mode = 'concat')
pm2= merge([cp3, cp4], mode= 'concat')
pm= merge([pm1, pm2], mode= 'concat')
pm= Dropout(0.1)(pm)
px = Flatten()(pm)

x = merge([sx, px], mode='concat')
x = Dense(64, init= 'he_uniform', activation='relu')(x)
x = Dense(64, activation='relu')(x)
x = Dropout(0.5)(x)

predictions = Dense(3, activation='softmax', name = '3_state')(x)
predictions_8 = Dense(8, activation = 'softmax', name = '8_state')(x)
predictions_rsa = Dense(1, activation = 'sigmoid', name = 'rsa')(x)

model_states = Model(input=[seq_input, pssm_input], output=[predictions, predictions_8])		
model_rsa = Model(input=[seq_input, pssm_input], output=predictions_rsa)		

model_states.compile(optimizer='adagrad',
              loss='categorical_crossentropy',
              metrics=['accuracy'])

model_rsa.compile(optimizer='adagrad',
              loss='mse',
              metrics=['accuracy'])
model_states.summary()
nb_epoch= 20
batchsize= 1000

history_states = model_states.fit([train_seq, train_pssm], [train_ss3, train_ss8], nb_epoch=nb_epoch, batch_size=batchsize, validation_data=([test_seq, test_pssm], [test_ss3, test_ss8]))
history_rsa = model_rsa.fit([train_seq, train_pssm], train_rsa, nb_epoch=nb_epoch, batch_size=batchsize, validation_data=([test_seq, test_pssm], test_rsa))

from keras.utils.visualize_util import plot
plot(model_states, to_file='model_states.png')
plot(model_rsa, to_file='model_rsa.png')

 

'''

loss=history.history['loss']    
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

#scores = model.evaluate(test_pssm, test_ss, batch_size=batchsize)
#print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100)), 'accuracy'    #prints accuracy, one value, in percent 

##############################################################################################################################
################################################### EVALUATION / PLOTTING ####################################################
##############################################################################################################################

#### confusion matrix
index = model.predict([test_ohe, test_pssm], batch_size= batchsize, verbose= 1)
index3= np.argmax(index[0], axis = 1)       # axis 1 gives you the max per row not total
index8= np.argmax(index[1], axis = 1)
true_val = test_ss.nonzero()[1]      # nonzero returns two items, we want the second
true_val_8 = test_ss_8.nonzero()[1]

### index and trueval should have the same numbers corresponding to each index. Helix is 0, sheet is 1, and coil is 2
## for 8 state: 'H':1, 'I':2 , 'G':3, 'E':4, 'B':5, 'T':6, 'S':7, '-':8
###H = a-helix, I = 5 helix (pi-helix) , G = 3-helix (310 helix), E = extended strand, participates in B ladder, B = residue in isolated B-bridge
### T = hydrogen bonded turn, S = bend   - = coil (residues which are not in any of the above conformations).


matrix = confusion_matrix(true_val, index3)
matrix_8 = confusion_matrix(true_val_8, index8)
print matrix
print matrix_8	
plt.figure()
plt.imshow(matrix)
plt.colorbar()
plt.figure()
plt.imshow(matrix_8)
plt.colorbar()
plt.show()


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

#### plotting predictions ### 
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

		
		print len(sheet),len(s_pos)
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
		
'''




### confusion matrix
		