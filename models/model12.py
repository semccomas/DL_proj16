import sys 
from keras.models import Model
from keras.layers import Dense, Dropout, Flatten, Convolution1D, MaxPooling1D, Input, merge
import numpy as np
import tables as tb
from sklearn.metrics import confusion_matrix, accuracy_score, precision_score, recall_score, classification_report 
import matplotlib.pyplot as plt
from keras.regularizers import l2, activity_l2
from mpl_toolkits.mplot3d import Axes3D


##############################################################################################################################
###################################################### LOADING TABLES ############################################################
##############################################################################################################################
def inputs (filename, seq, pssm, ss3, ss8, rsa, angles):
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
					angles.append(array.angles_feat.read()[where])
				except AttributeError:
					pass
	dataset.close()

######################## sequence info ##########################
table= 'big_table_all3'			
train_seq = []
train_pssm = []
train_ss3 = []
train_ss8 = []
train_rsa = []
train_angles = []
inputs(table, train_seq, train_pssm, train_ss3, train_ss8, train_rsa, train_angles)
train_seq= np.concatenate(train_seq, axis = 0)
train_pssm= np.concatenate(train_pssm, axis = 0)			
train_ss3= np.concatenate(train_ss3, axis = 0)			
train_ss8= np.concatenate(train_ss8, axis = 0)			
train_rsa= np.concatenate(train_rsa, axis = 0)			
train_angles = np.concatenate(train_angles, axis = 0)
	
test_table= 'big_test_table_all3'						
test_seq = []
test_pssm = []
test_ss3 = []
test_ss8 = []
test_rsa = []
test_angles = []
inputs(test_table, test_seq, test_pssm, test_ss3, test_ss8, test_rsa, test_angles)
test_seq= np.concatenate(test_seq, axis = 0)
test_pssm= np.concatenate(test_pssm, axis = 0)			
test_ss3= np.concatenate(test_ss3, axis = 0)			
test_ss8= np.concatenate(test_ss8, axis = 0)			
test_rsa= np.concatenate(test_rsa, axis = 0)	
test_angles = np.concatenate(test_angles, axis = 0)


def scale (angle_array):
	array = []
	for line in angle_array:
		line = line/360
		tup= (line[0], line[1])
		array.append(tup)
	array = np.asarray(array)
	return array 
train_angles= scale(train_angles)
test_angles = scale(test_angles)


train_rsa= np.nan_to_num(train_rsa)

print np.shape(test_angles), np.shape(train_angles), 'shape angles'
##############################################################################################################################
###################################################### THE MODEL ############################################################
##############################################################################################################################
seq_input = Input(shape=(20,15))
cs = Convolution1D(10, 5, border_mode = 'same',  activation= 'relu')(seq_input)
cs = Dropout(0.4)(cs)
cs = Convolution1D(10, 5,  activation= 'relu')(cs)
cs = Dropout(0.4)(cs)
#cs = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cs)
#cs = Dropout(0.4)(cs)
#cs = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cs)
#cs = Dropout(0.4)(cs)
#cs = Convolution1D(50, 5, border_mode = 'same', activation= 'relu')(cs)
#cs = Dropout(0.4)(cs)
#cs = Convolution1D(50, 5, border_mode = 'same', activation= 'relu')(cs)
#cs = Dropout(0.4)(cs)
cs = Convolution1D(10, 5, activation= 'relu')(cs)
cs = Dropout(0.4)(cs)
#cs = Convolution1D(50, 5, border_mode = 'same', activation= 'relu')(cs)
#cs = Dropout(0.4)(cs)
cs = Convolution1D(10, 5, activation= 'relu')(cs)
cs = Dropout(0.4)(cs)
cs = Convolution1D(10, 5, activation= 'relu')(cs)
sx = Flatten()(cs)

pssm_input = Input(shape=(21,15)) 
cp = Convolution1D(10, 5, border_mode = 'same', activation= 'relu')(pssm_input) 
cp = Dropout(0.4)(cp)
cp = Convolution1D(10, 5, activation= 'relu')(cp)
cp = Dropout(0.4)(cp)
#cp = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cp)
#cp = Dropout(0.4)(cp)
#cp = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cp)
#cp = Dropout(0.4)(cp)
#cp = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cp)
#cp = Dropout(0.4)(cp)
#cp = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cp)
#cp = Dropout(0.4)(cp)
cp = Convolution1D(10, 5, activation= 'relu')(cp)
cp = Dropout(0.4)(cp)
#cp = Convolution1D(50, 5, border_mode = 'same',  activation= 'relu')(cp)
#cp = Dropout(0.4)(cp)
cp = Convolution1D(10, 5, activation= 'relu')(cp)
cp = Dropout(0.4)(cp)
cp = Convolution1D(10, 5, activation= 'relu')(cp)
px = Flatten()(cp)

x = merge([sx, px], mode='concat')
x = Dropout(0.5)(x)
x = Dense(500, activation='relu', W_regularizer=l2(10 **(-6)))(x)
x = Dropout(0.5)(x)

predictions = Dense(3, init='he_uniform', activation='softmax', name = '3_state')(x)
predictions_8 = Dense(8, activation = 'softmax', name = '8_state')(x)
predictions_rsa = Dense(1, name = 'rsa')(x)		## == main aux output
predictions_angles = Dense(2, name = 'angles')(x)

model = Model(input=[seq_input, pssm_input], output=[predictions, predictions_8, predictions_rsa, predictions_angles])		

model.compile(optimizer='adam',  ## maybe change to adam
              loss={'3_state' :'categorical_crossentropy', '8_state' : 'categorical_crossentropy', 'rsa': 'mse', 'angles':'mse'} ,
              metrics={'3_state': 'accuracy', '8_state': 'accuracy'},
              loss_weights={'3_state' :1, '8_state' : 0.2, 'rsa': 1, 'angles': 1}) ##mse??
model.summary()
nb_epoch= 500
batchsize= 1000
from keras.utils.visualize_util import plot
plot(model, to_file='model.png', show_shapes=True)

history = model.fit([train_seq, train_pssm], {'3_state': train_ss3, '8_state': train_ss8, 'rsa': train_rsa, 'angles': train_angles}, nb_epoch=nb_epoch, batch_size=batchsize, verbose = 2, validation_data=([test_seq, test_pssm], {'3_state': test_ss3, '8_state': test_ss8, 'rsa': test_rsa, 'angles': test_angles}))
#exit()
 
loss=history.history['loss']    
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

#scores = model.evaluate(test_pssm, test_ss, batch_size=batchsize)
#print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100)), 'accuracy'    #prints accuracy, one value, in percent 


##############################################################################################################################
################################################### EVALUATION / PLOTTING ####################################################
##############################################################################################################################



############################################### confusion matrix#############################################################



pred = model.predict([test_seq, test_pssm], batch_size= batchsize, verbose= 1)
pred3= np.argmax(pred[0], axis = 1)       # axis 1 gives you the max per row not total
pred8= np.argmax(pred[1], axis = 1)
predR= (pred[2])
true_val3 = np.argmax(test_ss3, axis = 1)      # nonzero returns two items, we want the second
true_val8 = np.argmax(test_ss8, axis = 1)

### index and trueval should have the same numbers corresponding to each index. Helix is 0, sheet is 1, and coil is 2
## for 8 state: 'H':1, 'I':2 , 'G':3, 'E':4, 'B':5, 'T':6, 'S':7, '-':8
###H = a-helix, I = 5 helix (pi-helix) , G = 3-helix (310 helix), E = extended strand, participates in B ladder, B = residue in isolated B-bridge
### T = hydrogen bonded turn, S = bend   - = coil (residues which are not in any of the above conformations).

print 'CONFUSION MATRICES'
print
matrix = confusion_matrix(true_val3, pred3)
print matrix
matrix_8 = confusion_matrix(true_val8, pred8)
print matrix_8	
plt.figure()
plt.imshow(matrix)
plt.colorbar()
plt.figure()
plt.savefig('matrix3_model7')
plt.imshow(matrix_8)
plt.colorbar()
plt.savefig('matrix8_model7')
#plt.show()



############################################### classification report (precision and recall) #############################################################
#print np.cumsum(test_ss3, axis = 0)[-1], 'sum of each occurence, helix, sheet, coil (respectively)'
#print np.cumsum(test_ss8, axis = 0)[-1], 'sum of each occurence in 8 states, alpha helix[0], 5helix[1], 3helix[2], extended strand [3], b bridge[4], h-bonded turn[5], bend[6], and coil[7]'
print 
print 
print 
print 'CLASSIFICATION REPORT'
print '3 states'
print classification_report(true_val3, pred3)
print 
print '8 states'
print classification_report(true_val8, pred8)
print '(support ==  total number of occurences in each class) Nothing in class 1'
print 
print 
print 




############################################### triangle plot #############################################################


helix_mask = true_val3 == 0
fig = plt.figure()
ax= fig.add_subplot(111, projection = '3d')
x = pred[0][helix_mask,0]
y= pred[0][helix_mask,1]
z=pred[0][helix_mask,2]
ax.scatter(x,y,z, c = 'r', marker= 'o', alpha=0.1)
ax.scatter(1,0,0, c='k')
ax.scatter(0,1,0, c='k')
ax.scatter(0,0,1, c='k')
ax.plot([1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], c='k')
ax.set_xlabel('helix')
ax.set_ylabel('sheet')
ax.set_zlabel('coil')
#plt.savefig('helix_model7_3d')
plt.show()

sheet_mask = true_val3 == 1
fig = plt.figure()
ax= fig.add_subplot(111, projection = '3d')
x = pred[0][sheet_mask,0]
y= pred[0][sheet_mask,1]
z=pred[0][sheet_mask,2]
ax.scatter(x,y,z, c = 'g', marker= 'o', alpha=0.1)
ax.scatter(1,0,0, c='k')
ax.scatter(0,1,0, c='k')
ax.scatter(0,0,1, c='k')
ax.plot([1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], c='k')
ax.set_xlabel('helix')
ax.set_ylabel('sheet')
ax.set_zlabel('coil')
plt.show()
#plt.savefig('sheet_model7_3d')



coil_mask = true_val3 == 2
fig = plt.figure()
ax= fig.add_subplot(111, projection = '3d')
x = pred[0][coil_mask,0]
y= pred[0][coil_mask,1]
z=pred[0][coil_mask,2]
ax.scatter(x,y,z, c = 'b', marker= 'o', alpha=0.1)
ax.scatter(1,0,0, c='k')
ax.scatter(0,1,0, c='k')
ax.scatter(0,0,1, c='k')
ax.plot([1, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0], c='k')
ax.set_xlabel('helix')
ax.set_ylabel('sheet')
ax.set_zlabel('coil')
#plt.savefig('coil_model7_3d')
plt.show()



############################################### plotting loss #############################################################
epoch= []
for r in xrange(nb_epoch):
    epoch.append(r)

plt.figure(3)
plt.plot(epoch, loss)
plt.plot(epoch, history.history['val_loss'])
plt.ylabel('loss')
plt.xlabel('epoch')
plt.ylim(ymin=0)
#plt.savefig('loss_model7')
plt.show()

#####################################################################################
############################ measuring per protein ###################################
#####################################################################################
o = open('suckers', 'a')
alignment_lens= []
alignment_comp = open('line_len_test_aln').read().splitlines()

data= tb.open_file(test_table)
for group in data.walk_groups():
	for array in group:
		try:
			pssm = array.pssm_tab.read()
			seq = array.seq_tab.read()
			states3 = array.ss3_feat.read()
			states8 = array.ss8_feat.read()
			rsa = array.rsa_feat.read()

		except AttributeError:
			pass

		name= array._v_name
		

		predict=model.predict([seq, pssm], batch_size= batchsize, verbose= 1)
		maximum= np.amax(predict[0], axis= 1)				# max probability in each row
		positions3= np.argmax(predict[0], axis= 1)			# which index that (^) probability is. Is same as model.predict_classes
		positions8 = np.argmax(predict[1], axis = 1)
		actual_positions3 = np.argmax(states3, axis = 1)
		actual_positions8 = np.argmax(states8, axis = 1)
		
		print
		a3 = accuracy_score(actual_positions3, positions3)
		a8 = accuracy_score(actual_positions8, positions8)
		if a3 <= 0.5:
			print a3, 'accuracy score for 3 state for protein %r'  %(name)
			o.write(str(a3) + ' acc_3 ' + name + '\n')
			if a3 <= 0.2:
				print 'THIS ONE IS REALLY BAD!!!!!!!!!!!'
		if a8 <= 0.5:
			print a8, 'accuracy score for 8 states for protein %r' %(name)
			o.write(str(a8)+ ' acc_8 ' + name + '\n' + '\n')
			if a8 <= 0.2:
				print 'THIS ONE IS REALLY BAD!!!!!!!!!!!!!!'
		print
		

		### plotting acc vs aln len 

		for fline in alignment_comp:
			fline = fline.split()
			fname= fline[1][-14:-9]
			if name[-5:] in fname:
				tup = (a3, fline[0])
		alignment_lens.append(tup)




		# compares two lists, list 2 is an index for the 3 values in list one, appends the maximum value to the appropriate group. Rest is 0. 
		# pos is to later use in graph to show if predicted to be present. Really the only diff is threshold
		def find (lst, number, positions, actual):
			for index in xrange(len(positions)):
				if positions[index] == number:
					if actual == 1:
						lst.append(0.5)
					else:
						maximum= np.amax(predict[0], axis= 1)
						if maximum[index] >= 0.6:		#### !!!!!!!!!!!!! This is the one you can change if you want it to be a different threshold
							lst.append(maximum[index])		
						else:
							lst.append(0)
				else: 
					if actual == 1:
						lst.append(np.nan)   #nans arent plotted
					else:
						lst.append(0)


		helix= []
		h_pos= []
		sheet= []
		s_pos= []
		coil= []
		c_pos=[]		
		find(helix, 0, positions3, 0)
		find(sheet, 1, positions3, 0)
		find(coil, 2, positions3, 0)
		find(h_pos, 0, actual_positions3, 1)
		find(s_pos, 1, actual_positions3, 1)
		find(c_pos, 2, actual_positions3, 1)
		'''
	
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
		
		#plt.show()


		rsa_true = rsa
		rsa_true = np.nan_to_num(rsa_true)		 
		rsa_pred = predict[2]
		index = np.arange(len(rsa_true))
		bar_width = 0.35
		opacity = 0.4

		plt.bar(index, rsa_true, bar_width,
                 alpha=opacity,
                 color='b',
                 label='rsa_true')

		plt.bar(index + bar_width, rsa_pred, bar_width,
                 alpha=opacity,
                 color='r',
                 label='pred')

		print name 

		repeat = raw_input('Continue? [y/n]')
		print repeat
		if repeat == 'y':
			continue
		else:
			sys.exit()

		'''


o.close()


alignment_lens = np.asarray(alignment_lens)
a = sorted(alignment_lens, key=lambda alignment_lens_entry: len(alignment_lens_entry[1]))
a = np.asarray(a)
plt.plot(a[:,1], a[:,0])
plt.show()
