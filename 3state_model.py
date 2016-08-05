from keras.models import Sequential
from keras.layers import Dense, Dropout
#from keras.callbacks import TensorBoard
import numpy as np
import tables as tb

np.random.seed(3)
##############################################################################################################################
###################################################### LOADING TABLES ############################################################
##############################################################################################################################

############################### TRAINING SET ############################
train_ss = []
train_ohe = []
dataset = tb.open_file('big_table')
for group in dataset.walk_groups():
	for array in group:
		try:
			where = np.where(np.max(array.ss.read(), axis=1) == 1)
			#print where
			#raw_input()
			train_ss.append(array.ss.read()[where]) 
			train_ohe.append(array.one_hot.read()[where])
		except AttributeError:
			pass

train_ss = np.concatenate(train_ss, axis=0)
train_ohe = np.concatenate(train_ohe, axis = 0)

####################### TESTING SET #########################
test_ss = []
test_ohe = []
testset = tb.open_file('big_test_table')
for group in testset.walk_groups():
	for array in group:
		try:
			where = np.where(np.max(array.ss.read(), axis=1) == 1)
			#print where
			#raw_input()
			test_ss.append(array.ss.read()[where]) 
			test_ohe.append(array.one_hot.read()[where])
		except AttributeError:
			pass

test_ss = np.concatenate(test_ss, axis=0)
test_ohe = np.concatenate(test_ohe, axis = 0)

print len(train_ohe), len(train_ss)
print len(test_ohe), len(test_ss)

## change this later 
X_train= train_ohe
Y_train= train_ss
X_test= test_ohe
Y_test= test_ss


##############################################################################################################################
###################################################### THE MODEL ############################################################
##############################################################################################################################


model = Sequential()
model.add(Dense(64, input_dim=1, init='uniform', activation='tanh'))
#model.add(Dropout(0.5))
model.add(Dense(64, activation='tanh'))
#model.add(Dropout(0.5))
model.add(Dense(3, activation='softmax'))



model.compile(loss='categorical_crossentropy',
              optimizer='adagrad',
              metrics=['accuracy'])




## sparse categorical was raising errors so I changed to just categorical

batchsize= 1000         #using batch size more than once so defining it as a var to not forget to change all values each time

nb_epoch= 10
#TB= TensorBoard(log_dir='./logs', histogram_freq=0, write_graph=True)
history = model.fit(X_train, Y_train, nb_epoch=nb_epoch, batch_size=batchsize, validation_data=(X_test, Y_test))

print history.history
loss=history.history['loss']
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

scores = model.evaluate(X_test, Y_test, batch_size=batchsize)
print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))


#from keras.utils.visualize_util import plot
#plot(model, to_file='model.png')

##############################################################################################################################
###################################################### EVALUATION ############################################################
##############################################################################################################################

import matplotlib.pyplot as plt

x= []
for r in xrange(nb_epoch):
    x.append(r)
print x

plt.plot(x, loss)
plt.ylabel('loss')
plt.xlabel('epoch')
plt.show()
plt.savefig('old_dssp_loss.png')





dataset.close()
testset.close()

