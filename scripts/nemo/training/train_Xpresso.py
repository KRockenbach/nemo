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
title: train_Xpresso.py
description: trains Xpresso model for comparison with nemo
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-20
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/training/train_Xpresso.py <number of replications> <TPM type ["max" or "median"]>

notes: use nemo environment to run this script
=========================================================================================================
'''

import sys, h5py, os
import numpy as np
from math import ceil
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers.experimental import Adam
from tensorflow.keras.callbacks import EarlyStopping, CSVLogger, TerminateOnNaN, ModelCheckpoint, LambdaCallback, ReduceLROnPlateau
from tensorflow.keras.utils import plot_model
from tensorflow.keras.metrics import R2Score
from ..utils import model_utils
import tensorflow.keras.backend as K
from yaml import dump
from ..models.Xpresso import build_xpresso

print("\n\n")

root_dir = os.path.join('..') # relative to bash script

organism="Bnapus"
model_descriptor = "xpresso"
modeldir = os.path.join(root_dir, 'model_configs', model_descriptor)
resultdir = modeldir.replace('model_configs', 'model_weights')
orgdir = os.path.join(resultdir, organism)
mask_type = "masked"
partition_type = "graphpart_Bn"
outdir = os.path.join(orgdir, f'{mask_type}_{partition_type}')
datadir=os.path.join(root_dir, "data", organism, f'{mask_type}_{partition_type}_fold_data')

logging = True

test_fold = 0
valid_fold = 1

N = int(sys.argv[1]) # rep number

tpm_type=argv[2] # output: max or median


if tmp_type == "max" or tpm_type == "maximum":
    out_idx = 4
else:
    out_idx = 2

# create output directory
os.makedirs(outdir, exist_ok=True)

# define hyperparameter dictionary
hp_path = os.path.join(modeldir, "hyperparams.tsv") #tsv file containing hyperparameters
params = model_utils.dict_from_tsv(hp_path)

# cofig contains hyperparameters for different model types
# hyperparams in config are more general and not optimized during hyperparameter tuning
conf_path = os.path.join(modeldir, "config.tsv")
config = model_utils.dict_from_tsv(conf_path)

input_names = ["promoter", "halflife"]
outP = int(params["outside_P"])
inP = int(params["inside_P"])
outT = 0
inT = 0

# get data for the training and validation set depending on config
#########################

train = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
train["output"] = train["output"][,out_idx]
valid = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)
valid["output"] = valid["output"][,out_idx]
n_training = train[input_names[0]].shape[0]
print(f"{n_training} training examples")
print("\n\n")
########################


K.clear_session()
#build model
#model = model_utils.build_model(params=params, train=train, input_names=input_names)
model = build_xpresso()
print(model.summary())
print("\n\n")

if N == 1:
    #plot model
    plot_model(model,
            to_file=(modeldir + '/' + model_descriptor + '_model.png'),
            show_shapes=True,
            show_dtype=False,
            show_layer_names=True,
            rankdir='TB',
            expand_nested=False,
            dpi=96,
            layer_range=None,
            show_layer_activations=True,
            show_trainable=False
            )

    model_json = model.get_config()
    yaml.dump(model_json, os.path.join(modeldir, "keras_config.yaml"), allow_unicode=True)


#### callbacks ####
logdir = os.path.join(outdir, 'logs')
logdir = logdir.replace("model_weights", "results")
os.makedirs(logdir, exist_ok=True)

# val_loss needed for model checkpoints
model_outfile = f"{model_descriptor}_rep{N}_{tpm_type}.h5"
check_cb = ModelCheckpoint(os.path.join(outdir, model_outfile),
                           monitor='val_loss', verbose=1,
                           save_best_only=True, mode='min')

# stop if val_loss has not decresed for 10 epochs
earlystop_cb = EarlyStopping(monitor='val_loss', patience=10, verbose=1, mode='min')
logpath = os.path.join(logdir, f'trainlog_rep{N}_{tpm_type}.csv')


# terminate on NaN training loss
term_cb = TerminateOnNaN()

# Reduce learning rate on plateau
redLR_cb = ReduceLROnPlateau()

# redLR_cb, check_cb and earlystop_cb added later, depending on config (not needed when using one-cycle LR policy, or training on the full set)
# log metrics/losses at each epoch
csvlog_cb = CSVLogger(logpath, append=True, separator='\t')
cb_list = [csvlog_cb, term_cb]


#### compile model ####
optimizer_fxn = Adam(learning_rate=0.001,
                     beta_1=0.9,
                     beta_2=0.999,
                     epsilon=1e-08,
                     weight_decay=0.0)
loss_fxn = "mse"
print("Using Adam with MSE loss")

model.compile(optimizer_fxn,
              loss = loss_fxn,
              metrics=[R2Score()])

cb_list.append(check_cb)
cb_list.append(redLR_cb)
cb_list.append(earlystop_cb)
num_epochs=500 # maximum number of epochs

# train model
batch=int(float(params['batch_size']))

print(f"Training model with fold configuration test {str(test_fold)} valid {str(valid_fold)}")
print(f"Organism: {organism}")
model.fit([train[i] for i in input_names],
          train["output"],
          batch_size=batch,
          epochs=num_epochs,
          validation_data=([valid[i] for i in input_names],
                           valid["output"]),
          callbacks=cb_list,
          verbose=1) # display progress bars
del model





