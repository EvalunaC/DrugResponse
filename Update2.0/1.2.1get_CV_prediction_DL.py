# -*- coding: utf-8 -*-

"""
Created on Fri Nov  5 22:20:04 2021

@author: EvalunaC
"""

import pyreadr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf
import tensorflow_addons as tfa

#from sklearn.model_selection import cross_val_score
from sklearn.metrics import r2_score


from sklearn.model_selection import train_test_split
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.optimizers import Adam


from keras.optimizers import schedules
from keras.optimizers import SGD

print(tf.__version__)
print(np.__version__)
import platform
print(platform.python_version())

"""### Importing the training set"""

with open('/extraspace1/qli12/Data Source/GDSC2/192Drug_list.csv') as f:
    drugnames = [line.rstrip().strip('\"') for line in f]

train_file = "/extraspace1/qli12/Data Source/GDSC2/GDSC2_Frame_NoTrans/"+drugnames+"_trainFrame.RData"


    
    
def DL(drug):
    print(i)
    print("=============================")
    train_file = "/extraspace1/qli12/Data Source/GDSC2/GDSC2_Frame_NoTrans/"+drug+"_trainFrame.RData"
    trainFrame =pyreadr.read_r(train_file)['trainFrame']
    trainFrame.shape
    X = trainFrame.iloc[:,1:].values
    y = trainFrame.iloc[:,0].values

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = 0.2, random_state = 0)


#Variables
    layerOneNumNodes = 128
    layerTwoNumNodes = 128
    layerThreeNumNodes = 128

    activation = 'relu'

    learningRate = 0.01
    learningRateDecay = 1

    regressor = Sequential()
    regressor.add(Dense(units = layerOneNumNodes, input_shape = (X_train.shape[1],), activation = activation))
    regressor.add(Dropout(0.4))

    regressor.add(Dense(units = layerTwoNumNodes, activation = activation))
    regressor.add(Dropout(0.4))

    regressor.add(Dense(units = layerThreeNumNodes, activation = activation))
    regressor.add(Dropout(0.4))
    regressor.add(Dense(units = 1,  activation = 'linear'))

    lr_schedule = schedules.ExponentialDecay(
        initial_learning_rate=learningRate,
        decay_steps=10000,
        decay_rate=learningRateDecay)
    sgd = tf.keras.optimizers.SGD(
        learning_rate=learningRate, momentum=0.0, nesterov=False, name="SGD"
        )

    opt = Adam(learning_rate=lr_schedule)


    regressor.compile(optimizer = opt, loss = 'mse')

    history = regressor.fit(X_train, y_train, epochs = 50, batch_size = 16, shuffle=True)

# summarize regressor for loss
#    plt.plot(history.history['loss'][0:])
#    plt.title('model loss')
#    plt.ylabel('loss')
#    plt.xlabel('epoch')
#   plt.legend(['train', 'test'], loc='upper left')
#    plt.annotate(str(history.history['loss'][-1]), xy = (history.history['loss'][-1], 1),
#                 xytext = (.5, .5))

#save photo
    predictions = regressor.predict(X_test)
    predictions=np.reshape(predictions, (predictions.shape[0],))
    output = np.vstack((predictions,y_test,[drug]*predictions.shape[0])).T
    file = "/extraspace/ychen42/1.2.1.CV_perdiction_DL/"+drug+"_CV_DLpredict.csv"
    np.savetxt(file, output, fmt="%s", delimiter=',')


for i in range(1,len(drugnames)):
    DL(drugnames[i])
