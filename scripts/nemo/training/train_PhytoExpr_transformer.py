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
title: train_PhytoExpr_transformer.py
description: trains transformer model for comparison with nemo
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-20
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/training/train_PhytoExpr_transformer.py <number of replications> <TPM type ["max" or "median"]>

notes: use nemo environment to run this script
       based on work by Li et al. 2024 (https://doi.org/10.6084/m9.figshare.24417076, License: https://creativecommons.org/licenses/by/4.0/)
=========================================================================================================
'''

import numpy as np
import pandas as pd
from tensorflow.keras.models import Model
from tensorflow.keras.layers import *
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.optimizers import Adam
from sys import argv
import random
import os
import math
from ..models.PhytoExpr.transformer import *
from yaml import dump
from ..utils import model_utils
import tensorflow.keras.backend as K

N = int(sys.argv[1])

tpm_type=argv[2]             # median, max
test_fold=0
valid_fold=1
outside=4000
inside=1000
batch_size_train=64
batch_size_test=512

model_descriptor = "PhytoExpr_transformer"
modeldir = os.path.join(root_dir, 'model_configs', model_descriptor)
resultdir = modeldir.replace('model_configs', 'model_weights')
mask_type = "masked"
partition_type = "graphpart_Bn"
outdir = os.path.join(resultdir, "Bnapus", f'{mask_type}_{partition_type}')
datadir=os.path.join(root_dir, "data", organism, f'{mask_type}_{partition_type}_fold_data')


data_config = {'use_promoter': 'True', 'use_terminator': 'True', 'use_halflife': 'False'}

train_dict = model_utils.get_set(data_config, outP=outside, inP=inside, outT=outside, inT=inside, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
valid_dict = model_utils.get_set(data_config, outP=outside, inP=inside, outT=outside, inT=inside, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)

# concatenate promoter and terminator
train_input = np.concatenate(train_dict["promoter"], train_dict["terminator"], axis=1)
valid_input = np.concatenate(valid_dict["promoter"], valid_dict["terminator"], axis=1)

# get correct target output
if tpm_type == "max" or tpm_type == "maximum":
    out_idx = 4
else:
    out_idx = 2

train_output = train["output"] = train_dict["output"][,out_idx] # median or max
valid_output = train["output"] = valid_dict["output"][,out_idx] # median or max

del train_dict, valid_dict

K.clear_session()
###############################
# build model
input_length=2*(outside+inside)

input=Input(shape=(input_length, 4))
Procheck=procheck(input)
model=Model(inputs=input, outputs=Procheck)
model.summary()

total_params = model.count_params()

# training
callbacks=[EarlyStopping(monitor='val_loss',patience=2,verbose=0,restore_best_weights=True)]
logdir = os.path.join(outdir, 'logs')
logdir = logdir.replace("model_weights", "results")
os.makedirs(logdir, exist_ok=True)
logpath = os.path.join(logdir, f'trainlog_rep{N}.csv')
# log metrics/losses at each epoch
csvlog_cb = tf.keras,callbacks.CSVLogger(logpath, append=True, separator='\t')
callbacks.append(csvlog_cb)
adam=Adam(learning_rate=0.0001, beta_1=0.9, beta_2=0.999, decay=0.00, amsgrad=False)
model.compile(optimizer=adam,
              loss='mean_squared_error',
              metrics=['mse', R2Score()])
model.fit(train_input,
              train_output,
              batch_size=batch_size,
              epochs=100,
              validation_data=(valid_input,
                               valid_output),
              callbacks=callbacks,
              verbose=1) # display progress bars
model.save(os.path.join(outdir, f'TransformerModel_rep{N}_{tpm_type}.h5'))
if N==1:
    #plot model
    plot_model(model,
               to_file=(modeldir + '/transformer_model.png'),
               show_shapes=True,
               show_dtype=False,
               show_layer_names=True,
               rankdir='TB',
               expand_nested=False,
               dpi=200,
               layer_range=None,
               show_layer_activations=True,
               show_trainable=False
               )
    model_json = model.get_config()
    yaml.dump(model_json, os.path.join(modeldir, "keras_config.yaml"), allow_unicode=True)
del model



