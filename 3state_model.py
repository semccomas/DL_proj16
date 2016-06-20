## this is taken from the keraswork.py thing we started with. Just getting an idea of how we will make this first model

from keras.models import Sequential
from keras.layers import Dense, Dropout
from keras.callbacks import TensorBoard
import numpy as np


np.random.seed(3)

dataset = np.loadtxt("dssp_values", delimiter="\t")
X_train = dataset[:,0]
Y_train = dataset[:,1:]
testset=np.loadtxt('dssp_validation', delimiter='\t')
X_test = testset[:,0]
Y_test = testset[:,1:]



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

batchsize= 100         #using batch size more than once so defining it as a var to not forget to change all values each time


#TB= TensorBoard(log_dir='./logs', histogram_freq=0, write_graph=True)
history = model.fit(X_train, Y_train, nb_epoch=100, batch_size=batchsize, validation_data=(X_test, Y_test))

print history.history
#this will get us plottable numbers for making a graph. We have to include validation data into the history so that we actually have a testing set in there... 

scores = model.evaluate(X_test, Y_test, batch_size=batchsize)
print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))


#from keras.utils.visualize_util import plot
#plot(model, to_file='model.png')
