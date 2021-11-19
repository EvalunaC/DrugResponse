# -*- coding: utf-8 -*-

"""
Created on Fri Nov  5 22:20:04 2021

@author: EvalunaC
"""

"/extraspace/ychen42/validation-CTRDB/Epirubicin_trainFrame_GDSC2.csv"


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf
import tensorflow_addons as tfa

#from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score

print(tf.__version__)
print(np.__version__)
import platform
print(platform.python_version())

"""### Importing the training set"""
trainFrame = pd.read_csv("/extraspace/ychen42/validation-CTRDB/Epirubicin_trainFrame_GDSC2.csv")

testFrame = pd.read_csv("/extraspace/ychen42/validation-CTRDB/Epirubicin_testFrame_GDSC2.csv")

trainFrame.shape
testFrame.shape


X_train = trainFrame.iloc[:,2:].values
y_train = trainFrame.iloc[:,1].values

X_test = testFrame.iloc[:,1:].values


from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.optimizers import Adam

"""### Initialising the ANN"""


#Variables
layerOneNumNodes = 128
layerTwoNumNodes = 128
layerThreeNumNodes = 128

activation = 'relu'

learningRate = 0.01
learningRateDecay = 1

regressor = Sequential()



"""### Adding the first Dense layer and some Dropout regularisations"""

regressor.add(Dense(units = layerOneNumNodes, input_shape = (X_train.shape[1],), activation = activation))
regressor.add(Dropout(0.4))

"""### Adding a second Dense layer and some Dropout regularisation"""

regressor.add(Dense(units = layerTwoNumNodes, activation = activation))
regressor.add(Dropout(0.4))

"""### Adding a third Dense layer and some Dropout regularisation"""

regressor.add(Dense(units = layerThreeNumNodes, activation = activation))
regressor.add(Dropout(0.4))
"""### Adding the output layer"""

regressor.add(Dense(units = 1,  activation = 'linear'))



from keras.optimizers import schedules
from keras.optimizers import SGD
lr_schedule = schedules.ExponentialDecay(
    initial_learning_rate=learningRate,
    decay_steps=10000,
    decay_rate=learningRateDecay)
sgd = tf.keras.optimizers.SGD(
    learning_rate=learningRate, momentum=0.0, nesterov=False, name="SGD"
)

opt = Adam(learning_rate=lr_schedule)

"""### Compiling the RNN"""

regressor.compile(optimizer = opt, loss = 'mse')

print(regressor.summary())

"""### Fitting the RNN to the Training set"""

history = regressor.fit(X_train, y_train, epochs = 50, batch_size = 16, shuffle=True)


# summarize regressor for loss
plt.plot(history.history['loss'][0:])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.annotate(str(history.history['loss'][-1]), xy = (history.history['loss'][-1], 1),
             
            xytext = (.5, .5))

#save photo
predictions = regressor.predict(X_test)
np.savetxt('/extraspace/ychen42/validation-CTRDB/Epirubicin_testFrame_DLpredict.csv', predictions, delimiter=",")

