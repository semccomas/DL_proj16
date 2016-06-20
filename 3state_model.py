## this is taken from the keraswork.py thing we started with. Just getting an idea of how we will make this first model

from keras.models import Sequential
from keras.layers import Dense, Dropout
import numpy


numpy.random.seed(3)

dataset = numpy.loadtxt("dssp_values", delimiter=",")
X = dataset[:,0]
Y = dataset[:,1]

model = Sequential()
model.add(Dense(64, input_dim=1, init='uniform', activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(64, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(1, activation='sigmoid'))

model.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['accuracy'])


model.fit(X, Y, nb_epoch=10, batch_size=100)


scores = model.evaluate(X, Y)
print("%s: %.2f%%" % (model.metrics_names[1], scores[1]*100))


from keras.utils.visualize_util import plot
plot(model, to_file='model.png')
