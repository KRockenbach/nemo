#!/usr/bin/python3

##################################################################################
#
# MIT License
#
# Copyright (c) 2025 Kevin Rockenbach, Agnieszka Golicz
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
##################################################################################


'''
=========================================================================================================
title: train_PhytoExpr_CNN.py
description: trains ensemble CNN model for comparison with nemo
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-20
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/training/train_PhytoExpr_CNN.py <number of replications> <TPM type ["max" or "median"]>

notes: use nemo environment to run this script
       based on work by Li et al. 2024 (https://doi.org/10.6084/m9.figshare.24417076, License: https://creativecommons.org/licenses/by/4.0/)
=========================================================================================================
'''

import numpy as np
import pandas as pd
from keras.models import Model, load_model
from keras.layers import *
from keras.callbacks import EarlyStopping
from keras.utils import to_categorical
from sys import argv
import pickle
import random
import os
import math
from ..models.PhytoExpr.CNN_submodels import *
from yaml import dump

import keras.backend as K
from sklearn import linear_model
from sklearn.metrics import mean_squared_error
from ..utils import model_utils


N = int(sys.argv[1])

root=".."
tpm_type=argv[2]             # median, max
test_fold=0
valid_fold=1

model_descriptor = "PhytoExpr_CNN"
modeldir = os.path.join(root_dir, 'model_configs', model_descriptor)
resultdir = modeldir.replace('model_configs', 'model_weights')
mask_type = "masked"
partition_type = "graphpart_Bn"
outdir = os.path.join(resultdir, "Bnapus", f'{mask_type}_{partition_type}')
datadir=os.path.join(root_dir, "data", organism, f'{mask_type}_{partition_type}_fold_data')
os.makedirs(outdir, exist_ok=True)
os.makedirs(modeldir, exist_ok=True)

config = pd.read_csv(os.path.join(modeldir, "ensemble_model_cfg.csv"), delimiter=";", header='0')

model_types = config[:,"model_type"].tolist()
outside_intervals = config[:,"outside_interval"].tolist()
inside_intervals = config[:,"inside_interval"].tolist()
batch_size=128
param_df = config.iloc[:,4:-2]


K.clear_session()
x = []
submodel_dict={}
for submodel_index range(27):
    prefix = f"submodel{submodel_index}_"
    model_type = model_types[submodel_index] # C2D3, C3D3, C4D3, C2D2, C3D2, C4D2
    outside = outside_intervals[submodel_index]
    inside = inside_intervals[submodel_index]
    params = param_df.iloc[submodel_index,:].tolist()

    data_config = {'use_promoter': 'True', 'use_terminator': 'True', 'use_halflife': 'False'}

    train_dict = model_utils.get_set(data_config, outP=outside, inP=inside, outT=outside, inT=inside, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
    valid_dict = model_utils.get_set(data_config, outP=outside, inP=inside, outT=outside, inT=inside, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)

    # concatenate promoter and terminator
    train_input = np.concatenate(train_dict["promoter"], train_dict["terminator"], axis=1)
    valid_input = np.concatenate(valid_dict["promoter"], valid_dict["terminator"], axis=1)

    # get correct target output
    if tpm_type == "max" or tpm_type == "maximum":
        out_idx = 4 # maximum
    else:
        out_idx = 2 # median

    train_output = train_dict["output"][,out_idx]
    valid_output = valid_dict["output"][,out_idx]

    if submodel_index == 0:
        submodel_dict["ID"] = train_dict["ID"]
        submodel_dict["actual"] = train_output

    del train_dict, valid_dict

    # build the model
    input_length=2*(outside+inside)

    input=Input(shape=(input_length,4))
    exec(f"x{submodel_index}=build_{model_type}(input,parameters,prefix)")
    exec(f"x.append(x{submodel_index})")
    model=None
    exec(f"model=Model(inputs=input,outputs=x{submodel_index})")
    model.summary()

    total_params = model.count_params()
    flatten_width = model.get_layer('flatten').output_shape[1]

    # training
    callbacks=[EarlyStopping(monitor='val_loss',patience=2,verbose=0,restore_best_weights=True)]
    logdir = os.path.join(outdir, 'logs')
    logdir = logdir.replace("model_weights", "results")
    os.makedirs(logdir, exist_ok=True)
    logpath = os.path.join(logdir, f'trainlog_rep{N}.csv')
    # log metrics/losses at each epoch
    csvlog_cb = tf.keras,callbacks.CSVLogger(logpath, append=True, separator='\t')
    callbacks.append(csvlog_cb)

    model.compile(optimizer='adam',
                  loss='mean_squared_error',
                  metrics=['mse'])

    model.fit(train_input,
          train_output,
          batch_size=batch_size,
          epochs=100,
          validation_data=(valid_input,
                           valid_output),
          callbacks=callbacks,
          verbose=1) # display progress bars

    model.save(os.path.join(outdir, f'submodel{submodel_index}_rep{N}_{tpm_type}.h5'))
    submodel_dict[f"predictions{submodel_index}"] = model.predict(train_input, batch_size=batch_size)
    del model

x=Concatenate()(x)
x=Dense(1,name='second_layer_model_tpm')(x)

ensemble_model = Model(inputs=input,outputs=x)
submodels = []
for submodel_index in range(27):
    exec("submodel{submodel_index} = Model(inputs=ensemble_model.input,outputs=[ensemble_model.get_layer('submodel{submodel_index}_tpm_output').output])")
    exec("submodels.append(submodel{submodel_index})")


##########################################################################################################
# Weights for the second-layer model
##########################################################################################################

# submodel prediction data for the second-layer model
data = pd.DataFrame.from_dict(submodel_dict)

def weights_for_the_second_layer_model():
    x_train = data.iloc[:,2:]
    y_train = data.loc[:,'actual']

    model = linear_model.LinearRegression()
    model.fit(x_train,y_train)

    weights = [np.expand_dims(model.coef_,1),np.expand_dims(model.intercept_,1)]
    return weights

weights = weights_for_the_second_layer_model()
ensemble_model.get_layer('second_layer_model_tpm').set_weights(weights)
ensemble_model.get_layer('second_layer_model_tpm').trainable = False

##########################################################################################################
# Weights for the first layer models
##########################################################################################################

def get_weights_of_first_layer_model(path):
    model = load_model(path)
    weights = model.get_weights()
    del model
    return weights

# transfer weights
for model_index in range(27):
    print 'MODEL_INDEX: '+str(model_index)
    weights = get_weights_of_first_layer_model(os.path.join(outdir, f'submodel{submodel_index}_rep{N}_{tpm_type}.h5'))
    submodels[model_index].set_weights(weights)
    submodels[model_index].trainable = False

ensemble_model.save(os.path.join(outdir, f'EnsembleModel_rep{N}_{tpm_type}.h5'))
if N == 1:
    #plot model
    plot_model(model,
            to_file=(modeldir + '/ensemble_model.png'),
            show_shapes=True,
            show_dtype=False,
            show_layer_names=True,
            rankdir='TB',
            expand_nested=False,
            dpi=300,
            layer_range=None,
            show_layer_activations=True,
            show_trainable=False
            )
    model_json = ensemble_model.get_config()
    yaml.dump(model_json, os.path.join(modeldir, "keras_config.yaml"), allow_unicode=True)
del ensemble_model
