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
title: train_Basenji-5K.py
description: trains Basenji model for comparison with nemo
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-20
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/training/train_Basenji-5K.py <number of replications> <TPM type ["max" or "median"]>

notes: use nemo environment to run this script
       based on work by Qiu et al. 2025 (https://zenodo.org/records/14898182, License: MIT)
=========================================================================================================
'''

import sys, h5py, os
import numpy as np
from math import ceil
import pandas as pd
import tensorflow as tf
from tensorflow.keras.utils import plot_model
from ..utils import model_utils
from ..models.Basenji-5k import *
import tensorflow.keras.backend as K
from tensorflow.keras.metrics import R2Score


root_dir = '..' # relative to bash script

organism="Bnapus"
model_descriptor = "xpresso"
modeldir = os.path.join(root_dir, 'model_configs', model_descriptor)
resultdir = modeldir.replace('model_configs', 'model_weights')
orgdir = os.path.join(resultdir, organism)
mask_type = "masked"
partition_type = "graphpart_Bn"
outdir = os.path.join(orgdir, f'{mask_type}_{partition_type}')
datadir=os.path.join(root_dir, "data", organism, f'{mask_type}_{partition_type}_fold_data')

test_fold = 0
valid_fold = 1

N = int(sys.argv[1])

K.clear_session()
input_names = ["promoter"]
outP = 2500
inP = 2500
outT = 0
inT = 0

tpm_type = sys.argv[2]

config = {"use_promoter": 'True', "use_terminator": 'False', "use_halflife": 'False'}


# get data for the training and validation set depending on config
#########################
train = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
train["output"] = train["output"][,2] # median
valid = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)
valid["output"] = valid["output"][,2] # median
n_training = train[input_names[0]].shape[0]
print(f"{n_training} training examples")
print("\n\n")
########################

# create output directory
os.makedirs(outdir, exist_ok=True)
os.makedirs(modeldir, exist_ok=True)


def decay(epoch):
  if epoch < 3:
    return 1e-3

  else:
      return 1e-5


callbacks = [tf.keras.callbacks.ModelCheckpoint(filepath=os.path.join(outdir, f"Basenji-5K_rep{N}_{tpm_type}.h5"),
                                        monitor='val_loss',
                                        verbose=0,
                                        mode='auto' ,
                                        # save_weights_only=False,
                                        save_best_only= True),
            tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=5, verbose=0),
            tf.keras.callbacks.LearningRateScheduler(decay)]

logdir = os.path.join(outdir, 'logs')
logdir = logdir.replace("model_weights", "results")
os.makedirs(logdir, exist_ok=True)
logpath = os.path.join(logdir, f'trainlog_rep{N}_{tpm_tpye}.csv')
# log metrics/losses at each epoch
csvlog_cb = tf.keras,callbacks.CSVLogger(logpath, append=True, separator='\t')
callbacks.append(csvlog_cb)


K.clear_session()
model = basenji_model(input_shape=(5000,4), W=15,L=3)
model.compile(loss=tf.keras.losses.MeanSquaredError(),
              optimizer=tf.keras.optimizers.Adam(),
              metrics=R2Score())


batch_size = 32
epochs = 100

print(f"Training rep {N} of Basenji-5K, predicitng {tpm_type}")
print(f"Organism: {organism}")

if N == 1:
    #plot model
    plot_model(model,
            to_file=(modeldir + '/Basenji-5K_model.png'),
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


model.fit(train["promoter"],
          train["output"][:,5], # maximum expression
          batch_size=batch_size,
          epochs=epochs,
          validation_data=(valid["promoter"],
                           valid["output"][:,5]),
          callbacks=callbacks,
          verbose=1) # display progress bars
del model
