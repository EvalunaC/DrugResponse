# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 10:17:57 2021

@author: EYC
"""


from __future__ import print_function
import numpy as np
import tensorflow as tf
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn import preprocessing
from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor
from xgboost import XGBRegressor
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import tensorflow as tf
import tensorflow_addons as tfa


#file = sys.argv[1]
#forest_type = sys.argv[2]
dataset = pd.read_csv('/extraspace/ychen42/Drug_Response/Data/Qiangs_data_DL/Lapatinib_train_data_BRCA_NoTrans.csv')

expression = dataset.iloc[:,2:].values
response = dataset.iloc[:,1].values
gene_name = dataset.iloc[:,0].values

scaler = preprocessing.StandardScaler()
scaler.fit(expression)
expression = scaler.transform(expression)
n_genes = np.shape(expression)[1]

## train-test split
expression_train, expression_test, y_train, y_test = train_test_split(expression, response, test_size=0.2, shuffle=True)


n_trees = 1000


    xgb = XGBRegressor(n_estimators=n_trees,n_jobs=-1)
    
    xgb.fit(expression_train, y_train)

    f_importance = xgb.feature_importances_

    forest = xgb.get_booster().trees_to_dataframe()
    forest = list(forest.groupby(by='Tree'))
    forest = list(list(zip(*forest))[1])

    def recode_feature_name_array(f_nparray):
        def recode_feature_name(name_str):
            if name_str[0] is "f":
                return int(name_str[1:])
            else:
                return int(-2)
        return list(map(recode_feature_name, f_nparray))

    def parse_booster(tree_df):
        f_idx = recode_feature_name_array(np.array(tree_df['Feature']))
        nodes = np.array(tree_df['Node'])
        roots = np.array(tree_df.loc[tree_df['Feature']!='Leaf','Node'])
        right = np.array(nodes[2::2])
        left = np.array(nodes[1::2])
        edge_list = np.stack((roots, left, roots, right), axis=1)
        edge_list = np.reshape(np.take(f_idx, edge_list), [-1, 2])
        edge_list = edge_list[edge_list.min(axis=1)>=0,:]
        return edge_list


    elist = np.unique(np.vstack(list(map(parse_booster, forest))), axis=0)




    rf = RandomForestRegressor(n_estimators=n_trees, bootstrap=False, n_jobs=-1)
    rf.fit(expression_train, y_train)

    ## feature importance from the forest
    f_importance = rf.feature_importances_


    def parse_tree(decision_tree):
        tree = decision_tree.tree_
        parse_list = np.array(list(zip(tree.feature, tree.children_left, tree.children_right)))
        roots = np.array(range(np.shape(parse_list)[0]))
        edge_list = np.stack((roots, parse_list[:,1], roots, parse_list[:,2]), axis=1)
        edge_list = np.reshape(np.take(parse_list[:,0], edge_list), [-1, 2])
        edge_list = edge_list[edge_list.min(axis=1)>=0,:]
        return edge_list


    elist = np.unique(np.vstack(list(map(parse_tree, rf.estimators_))), axis=0)


selected = np.unique(elist)
x_train = expression_train[:, selected]
x_test = expression_test[:, selected]
x = expression[:, selected]
y = response

from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.optimizers import Adam
import pickle

"""### Initialising the pretrained encoder"""

#Variables
layerOneNumNodes = 128
layerTwoNumNodes = 128
layerThreeNumNodes = 128

activation = 'relu'

learningRate = 0.01
learningRateDecay = 1

regressor = Sequential()

regressor.add(Dense(units = layerOneNumNodes, input_shape = (x_train.shape[1],), activation = activation))
regressor.add(Dropout(0.3))

"""### Adding a second Dense layer and some Dropout regularisation"""

regressor.add(Dense(units = layerTwoNumNodes, activation = activation))
regressor.add(Dropout(0.4))

"""### Adding a third Dense layer and some Dropout regularisation"""

regressor.add(Dense(units = layerThreeNumNodes, activation = activation))
regressor.add(Dropout(0.4))
"""### Adding the output layer"""

regressor.add(Dense(units = 1,  activation = 'linear'))

"""Adding reduced learning rate"""

from keras.optimizers import schedules
from keras.optimizers import SGD
lr_schedule = schedules.ExponentialDecay(
    initial_learning_rate=learningRate,
    decay_steps=10000,
    decay_rate=learningRateDecay)
sgd = tf.keras.optimizers.SGD(
    learning_rate=learningRate, momentum=0.0, nesterov=False, name="SGD"
)

"""Adding the optomizer"""

opt = Adam(learning_rate=lr_schedule)

"""### Compiling the RNN"""

regressor.compile(optimizer = opt, loss = 'mse')

print(regressor.summary())

"""### Fitting the RNN to the Training set"""

history = regressor.fit(x, y, epochs = 50, batch_size = 16, shuffle=True)

result = regressor.evaluate(x_test, y_test)


dataset_test = pd.read_csv('/extraspace/ychen42/Drug_Response/Data/Qiangs_data_DL/Lapatinib_testFrame_BRCA_NoTrans.csv') # path to the TCGA that we want to get predictions from
# Delete ID column
names = dataset_test.iloc[:,0]
dataset_test = dataset_test.iloc[:,1:].values
dataset_test.shape
x_test = dataset_test[:, selected]

print(dataset_test.shape[0], dataset_test.shape[1])
#dataset_test = np.reshape(dataset_test, (1,dataset_test.shape[0], dataset_test.shape[1]))
predictions = regressor.predict(x_test)
y_hat = regressor.predict(x)

import numpy

numpy.savetxt('/extraspace/ychen42/Drug_Response/DeepLearning/ForgeNet_predict.csv', predictions, delimiter=",")


"""
Evaluation
"""

#RMSE
rmse = tf.keras.metrics.RootMeanSquaredError()
rmse.update_state(y, y_hat)
rmse.result().numpy()
# lapatinib 0.78914607

#MAE
mae = tf.keras.metrics.MeanAbsoluteError()
mae.update_state(y, y_hat, sample_weight=None)
mae.result().numpy()
# lapatinib  0.6305025

#R2
y_hat=np.reshape(y_hat, (y_hat.shape[0],))
R2 =  tfa.metrics.r_square.RSquare()
R2.update_state(y, y_hat)

R2.result().numpy()
# lapatinib 0.5605377

import pyreadr
selected = pd.read_csv('/extraspace/ychen42/Drug_Response/DeepLearning/Boruta_FNN/selected_gene.csv') # path to the TCGA that we want to get predictions from

data = pyreadr.read_r('/extraspace/ychen42/Drug_Response/Data/Qiangs_data_DL/Lapatinib_trainFrame_BRCA_NoTrans.RData') 
response = data.iloc[:,1].values


model = Sequential([
  Dense(64, activation='relu', input_shape=(selected.shape[1],)),
  Dense(64, activation='relu'),
  Dense(10, activation='softmax'),
])

model.compile(
  optimizer='adam',
  loss='mse'
#  metrics=['accuracy'],
)
model.fit(
  selected, # training data
  response, # training targets
  epochs=50,
  batch_size=32,
)












import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()

weights = {
    'h1': tf.Variable(tf.truncated_normal(shape=[n_features, n_hidden_1], stddev=weights_init_sd)),
    'h2': tf.Variable(tf.truncated_normal(shape=[n_hidden_1, n_hidden_2], stddev=weights_init_sd)),
    'h3': tf.Variable(tf.truncated_normal(shape=[n_hidden_2, n_hidden_3], stddev=weights_init_sd)),
    # 'h4': tf.Variable(tf.truncated_normal(shape=[n_hidden_3, n_hidden_4], stddev=weights_init_sd)),
    'out': tf.Variable(tf.truncated_normal(shape=[n_hidden_3, ], stddev=weights_init_sd))
}

biases = {
    'b1': tf.Variable(tf.constant(biases_init_value, shape=[n_hidden_1])),
    'b2': tf.Variable(tf.constant(biases_init_value, shape=[n_hidden_2])),
    'b3': tf.Variable(tf.constant(biases_init_value, shape=[n_hidden_3])),
    # 'b4': tf.Variable(tf.constant(biases_init_value, shape=[n_hidden_4])),
    'b_out': tf.Variable(tf.constant(biases_init_value, shape=[1]))
}



def fully_connected_layer(input, weight, bias, keep_prop, activation="relu"):
    layer = tf.add(tf.matmul(input, weight), bias)
    if activation is "tanh":
        layer = tf.nn.tanh(layer)
    elif activation is "sigmoid":
        layer = tf.nn.sigmoid(layer)
    else:
        layer = tf.nn.relu(layer)
    layer = tf.nn.dropout(layer, keep_prob=keep_prop)
    return layer


x_ = tf.placeholder(tf.float32, [None, n_features])
y_ = tf.placeholder(tf.float32, [None, 1])
keep_prob_ = tf.placeholder(tf.float32)
lr_ = tf.placeholder(tf.float32)
 




layer_1 = tf.add(tf.matmul(x_, tf.multiply(weights['h1'], partition)), biases['b1'])
layer_1 = tf.nn.relu(layer_1)
layer_2 = fully_connected_layer(layer_1, weights['h2'], biases['b2'], keep_prob_)
layer_3 = fully_connected_layer(layer_2, weights['h3'], biases['b3'], keep_prob_)
out = tf.add(tf.matmul(layer_3, weights['out']), biases['b_out'])

cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits_v2(logits=out, labels=y_))
optimizer = tf.train.AdamOptimizer(learning_rate=lr_).minimize(cost)
y_score = tf.nn.softmax(logits=out) # trainin prediction

## feature importance part of the computational graph
var_target = tf.reduce_sum(tf.abs(tf.multiply(weights['h1'], partition)), 0)
var_source = tf.reduce_sum(tf.abs(tf.multiply(weights['h1'], partition)), 1)
var_fully = tf.reduce_sum(tf.abs(weights['h2']), 1)
gedfn_importance = var_target + var_source + var_fully
nn_imp = tf.sparse_tensor_to_dense(tf.SparseTensor(np.reshape(selected, [-1, 1]), gedfn_importance, [len(f_importance)]))
nn_imp = tf.zeros([len(f_importance)]) + nn_imp

var_importance = nn_imp







## initiate training logs
loss_rec = np.zeros([training_epochs, 1])
training_eval = np.zeros([training_epochs, 1])
testing_eval = np.zeros([training_epochs, 1])

with tf.Session() as sess:

    sess.run(tf.global_variables_initializer())
    total_batch = int(np.shape(x_train)[0] / batch_size)

    ## Training cycle
    for epoch in range(training_epochs):

        avg_cost = 0.
        x_tmp, y_tmp = shuffle(x_train, y_train)
        # Loop over all batches
        for i in range(total_batch-1):

            batch_x, batch_y = x_tmp[i*batch_size:i*batch_size+batch_size], \
                               y_tmp[i*batch_size:i*batch_size+batch_size]
            _, c = sess.run([optimizer, cost], feed_dict={x_: batch_x, y_: batch_y,
                                                          keep_prob_: 1 - dropout,
                                                          lr_: learning_rate})
            # Compute average loss
            avg_cost += c / total_batch

        del x_tmp, y_tmp

        ## Display logs per epoch step
        if epoch % display_step == 0:

            ## Monitor training
            loss_rec[epoch] = avg_cost
            y_s = sess.run([y_score], feed_dict={x_: x_train, y_: y_train, keep_prob_: 1})
            y_s = np.reshape(np.array(y_s), [np.shape(x_train)[0], 2])[:, 1]
            acc = metrics.accuracy_score(y_train[:, 1], y_s > 0.5)
            auc = metrics.roc_auc_score(y_train[:, 1], y_s)
            training_eval[epoch] = [acc, auc]

            print("Epoch:", '%d' % (epoch), "cost =", "{:.9f}".format(avg_cost),
                  "Training accuracy:", round(acc,3), " Training auc:", round(auc,3))

    ## Testing
    y_s = sess.run([y_score], feed_dict={x_: x_test, y_: y_test, keep_prob_: 1})
    y_s = np.reshape(np.array(y_s), [np.shape(x_test)[0], 2])[:, 1]
    acc = metrics.accuracy_score(y_test[:, 1], y_s > 0.5)
    auc = metrics.roc_auc_score(y_test[:, 1], y_s)

    print("*****=====", "Testing accuracy: ", round(acc, 3), " Testing auc: ", round(auc, 3), "=====*****")

    var_imp = sess.run([var_importance])
    var_imp = np.reshape(var_imp, [len(f_importance)])

np.savetxt(file.split("/")[-1].split(".")[0]+"_forgeNetImportance.csv", var_imp, delimiter=",")


    for epoch in range(training_epochs):
            ## Monitor training
            loss_rec[epoch] = avg_cost
            y_s = sess.run([y_score], feed_dict={x_: x_train, y_: y_train, keep_prob_: 1})
            y_s = np.reshape(np.array(y_s), [np.shape(x_train)[0], 2])[:, 1]
#            acc = metrics.accuracy_score(y_train[:, 1], y_s > 0.5)
#            auc = metrics.roc_auc_score(y_train[:, 1], y_s)
#            training_eval[epoch] = [acc, auc]

            print("Epoch:", '%d' % (epoch), "cost =", "{:.9f}".format(avg_cost),
                  "Training accuracy:", round(acc,3), " Training auc:", round(auc,3))
