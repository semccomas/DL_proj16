import matplotlib.pyplot as plt
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
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
				#print where
				#raw_input()
				ohe.append(array.one_hot.read()[where])		### changed from array.one_hot.read() for pssm, one more down low to change too
				ss.append(array.ss.read()[where]) 
			except AttributeError:
				pass

#### this is what you change to test different things ####
table= '8state_table'
test_table= '8state_test_table'


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


 
X_train= train_ohe
Y_train= train_ss
X_test= test_ohe
Y_test= test_ss


print np.shape(Y_train)
print 

##############################################################################################################################
###################################################### THE MODEL ############################################################
##############################################################################################################################

'''
model = Sequential()
model.add(Convolution1D(64, 3, input_shape=(21, 15), input_length= 15, activation ='relu'))
model.add(Dropout(0.5))
model.add(Convolution1D(64, 3, activation= 'relu')) #, input_shape=(21, 15), input_length= 15 ))
model.add(Dropout(0.5))
model.add(Convolution1D(64, 3, activation= 'relu')) #, input_shape=(21, 15), input_length= 15 ))
model.add(Dropout(0.5))

'''

model = Sequential()
model.add(Flatten(input_shape=(20, 15)))		# This was changed from 20,15 to 21, 15
model.add(Dense(32, init='uniform', activation='tanh'))   
model.add(Dropout(0.5))
model.add(Dense(32, activation='tanh'))
model.add(Dropout(0.5))
model.add(Dense(8, activation='softmax'))



model.compile(loss='categorical_crossentropy',
              optimizer='adagrad',
              metrics=['accuracy'])



## sparse categorical was raising errors so I changed to just categorical

batchsize= 1000         #using batch size more than once so defining it as a var to not forget to change all values each time

nb_epoch= 20
#TB= TensorBoard(log_dir='./logs', histogram_freq=0, write_graph=True)
history = model.fit(X_train, Y_train, nb_epoch=nb_epoch, batch_size=batchsize, validation_data=(X_test, Y_test))

print history.history   #prints all values from the training periods (accuracy, val_accuracy, etc.)
loss=history.history['loss']    
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

scores = model.evaluate(X_test, Y_test, batch_size=batchsize)
print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100)), 'accuracy'    #prints accuracy, one value, in percent 



#predict=model.predict(X_test, batch_size= batchsize, verbose= 1)   #classify for the ss for each value. len == X_test len
#np.savetxt('bleepbloop', predict)



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


'''
dataset = tb.open_file(test_table)
for group in dataset.walk_groups():
	for array in group:
		single_protein= array.one_hot.read()		## this I also changed from array.one_hot.read() to pssm. Should only be two instances in whole file
		print np.shape(single_protein)
		predict=model.predict(single_protein, batch_size= batchsize, verbose= 1)
	#	raw_input()

		maximum= np.amax(predict, axis= 1)				# max probability in each row
		positions= np.argmax(predict, axis= 1)			# which index that (^) probability is. Is same as model.predict_classes
		raw_input()

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
		name= array._v_name
		print name

H = a-helix
B = residue in isolated B-bridge
E = extended strand, participates in B ladder
G = 3-helix (310 helix)
I = 5 helix (pi-helix)
T = hydrogen bonded turn
S = bend
'''