#this is from Keras- MLP for multi class softmax classification
#I think this is the best way to start for our 3 state model


from keras.models import Sequential
from keras.layers import Dense, Dropout, Activation
from keras.optimizers import SGD
import numpy as np


np.random.seed(3)
dataset= np.loadtxt('dssp_values', delimiter='\t')
X_train = dataset[:,0]
y_train = dataset[:,1:]


testset=np.loadtxt('dssp_validation', delimiter='\t')
X_test = testset[:,0]
y_test = testset[:,1:]




model = Sequential()
# Dense(64) is a fully-connected layer with 64 hidden units.
# in the first layer, you must specify the expected input data shape:
# here, 20-dimensional vectors.
model.add(Dense(64, input_dim=1, init='uniform'))
model.add(Activation('tanh'))
model.add(Dropout(0.5))
model.add(Dense(64, init='uniform'))
model.add(Activation('tanh'))
model.add(Dropout(0.5))
model.add(Dense(3, init='uniform'))
model.add(Activation('softmax'))

model.compile(loss='sparse_categorical_crossentropy',
              optimizer='adagrad',
              metrics=['accuracy'])

model.fit(X_train, y_train,
          nb_epoch=20,
          batch_size=160)
score = model.evaluate(X_test, y_test, batch_size=160)











########################## you want the output layer to be 3!! (Dense(3... because here we have 3 classes. We need to change the thing to be able to know that it is 3 classes not one

#maybe we have to change from softmax too 

#i changed model 1 input to 1, output to 1, and loss to be sparse_categorical instead of just categorical because keras said this expects a binary matrix only, no intergers with categorical. I can either convert or use sparse





#from keras.callbacks import Plotter
# standard plot
#model.fit(X_train, y_train, callbacks=[Plotter()])




