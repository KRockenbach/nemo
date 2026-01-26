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
title: run_training.py
description: trains model for a given species, masking-method and partitioning method, while holding out defined test/validation sets
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/run_training.py <species> <model_name> <masking> <partitioning> <test_fold_number> <valid_fold_number> <repetition_number>

notes: use nemo environment to run this script
       to run training on full set, set test_fold_number and valid_fold_number both to "None"
       to run training on 90% set, set only valid_fold_number to "None"
       logging and optimizer are set automatically depending on inputs
=========================================================================================================
'''

import sys, h5py, os
import numpy as np
from math import ceil
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.callbacks import EarlyStopping, CSVLogger, TerminateOnNaN, ModelCheckpoint, LambdaCallback, ReduceLROnPlateau
from tensorflow.keras.utils import plot_model
from tensorflow.keras.losses import Huber
from tensorflow.keras.metrics import R2Score
from ..utils import model_utils
import tensorflow.keras.backend as K

# important to use legacy version of Nadam for OneCycle callback to work
from tensorflow.keras.optimizers.legacy import Nadam
from ..utils.one_cycle_scheduler_tf.one_cycle_tf.one_cycle_scheduler import OneCycle
from tensorflow_addons.optimizers.weight_decay_optimizers import *
from yaml import dump
from ..models.nemo import build_nemo

print("\n\n")

root_dir = os.path.join('..') # relative to bash script

organism=sys.argv[1]
model_descriptor = "nemo"
modeldir = os.path.join(root_dir, 'model_configs', model_descriptor)
resultdir = modeldir.replace('model_configs', 'model_weights')
orgdir = os.path.join(resultdir, organism)
mask_type = sys.argv[2]
partition_type = sys.argv[3]
outdir = os.path.join(orgdir, f'{mask_type}_{partition_type}')
datadir=os.path.join(root_dir, "data", organism, f'{mask_type}_{partition_type}_fold_data')


logging = False
if partition_type == "graphpart_Bn":
    logging = True


if sys.argv[4] != "None":
    test_fold = int(sys.argv[4])
    full = False
    if sys.argv[5] != "None":
        valid_fold = int(sys.argv[5])
    else:
        valid_fold = None

else:
    full = True
    test_fold = None
    valid_fold = None
if sys.argv[6] == "None":
    N = None
else:
    N = int(sys.argv[6])

tpm_type = sys.argv[7]
if tpm_type == "max" or tpm_type == "maximum":
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



K.clear_session()
input_names = ["promoter", "terminator"]
outP = int(params["outside_P"])
inP = int(params["inside_P"])
outT = int(params["outside_T"])
inT = int(params["inside_T"])


if full:
    # get data for the training and validation set depending on config
    #########################
    train = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="full", test_fold=test_fold, valid_fold=valid_fold)
    train["output"] = train["output"][:,out_idx]
    n_training = train[input_names[0]].shape[0]
    print(f"{n_training} training examples")
    print("\n\n")
    ########################
else:
    # get data for the training and validation set depending on config
    #########################
    train = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
    train["output"] = train["output"][:,out_idx]
    valid = model_utils.get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)
    if valid_fold is not None:
        valid["output"] = valid["output"][:,out_idx]
    n_training = train[input_names[0]].shape[0]
    print(f"{n_training} training examples")
    print("\n\n")
    ########################



#build model
#model = model_utils.build_model(params=params, train=train, input_names=input_names)
#exec(f"model = model_utils.build_{model_descriptor}()")
model = build_nemo()
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
if logging:
    logdir = os.path.join(outdir, 'logs')
    logdir = logdir.replace("model_weights", "results")
    os.makedirs(logdir, exist_ok=True)
if full:
    model_outfile = f"{model_descriptor}_full.h5"
    if logging:
        # log metrics/losses at each epoch
        logpath = os.path.join(logdir, 'trainlog_full.csv')
    if N is not None:
        print(f"repeat: {N}")
        model_outfile = f"{model_descriptor}_full_n_{N}.h5"
        if logging:
            logpath = os.path.join(logdir, f'trainlog_full_n_{N}.csv')
    if logging:
        lrpath = logpath.replace("trainlog", "lr_log")
elif valid_fold is None:
    model_outfile = f"{model_descriptor}_t_{str(test_fold)}.h5"
    if N is not None:
        print(f"repeat: {N}")
        model_outfile = model_outfile.replace(".h5", f"_n_{N}.h5")
    if logging:
        # log metrics/losses at each epoch
        logpath = os.path.join(logdir, f"trainlog_t_{str(test_fold)}.csv")
        lrpath = logpath.replace("trainlog", "lr_log")
else:
    # val_loss needed for model checkpoints
    model_outfile = f"{model_descriptor}_t_{str(test_fold)}_v_{str(valid_fold)}.h5"
    if N is not None:
        print(f"repeat: {N}")
        model_outfile = model_outfile.replace(".h5", f"_n_{N}.h5")
    check_cb = ModelCheckpoint(os.path.join(outdir, model_outfile),
                           monitor='val_loss', verbose=1,
                           save_best_only=True, mode='min')

    # stop if val_loss has not decresed for 10 epochs
    earlystop_cb = EarlyStopping(monitor='val_loss', patience=10, verbose=1, mode='min')
    if logging:
        logpath = os.path.join(logdir, 'trainlog_t_{}_v_{}.csv'.format(str(test_fold), str(valid_fold)))
        lrpath = logpath.replace("trainlog", "lr_log")


# terminate on NaN training loss
term_cb = TerminateOnNaN()

# Reduce learning rate on plateau
redLR_cb = ReduceLROnPlateau()

# redLR_cb, check_cb and earlystop_cb added later, depending on config (not needed when using one-cycle LR policy, or training on the full set)
if logging:
    # log metrics/losses at each epoch
    csvlog_cb = CSVLogger(logpath, append=True, separator='\t')
    cb_list = [csvlog_cb, term_cb]
else:
    cb_list = [term_cb]



#### compile model ####

# one_cycle optimization

#max_lr  0.0663012475147729
#Bnorm_momentum  0.816691400963907
#huber_delta     2.9465947350775075
cycle_epochs = int(params['epochs'])
shift_peak = float(params['shift_peak'])
final_lr_scale = float(params['final_lr_scale'])
batch = int(params['batch_size'])
steps_per_epoch = ceil(n_training/float(batch))
cycle_size = cycle_epochs * steps_per_epoch
if logging:
    lrpath = lrpath.replace("lr_log", f"lr_log_{steps_per_epoch}_{n_training}")

max_lr = float(params['max_lr'])
min_lr = float(params['min_lr'])

# Nadam_beta_1 set via scheduler
Nadam_1_minus_beta_2 = float(params['Nadam_1_minus_beta_2'])
Nadam_beta_2 = 1 - Nadam_1_minus_beta_2
#lr_over_epsilon = 0.04872829743137629 / 8.327220255431336e-07
Nadam_epsilon = float(params['Nadam_epsilon'])
Nadam_weight_decay = float(params['Nadam_weight_decay'])

LR_schedule = OneCycle(
            initial_learning_rate=min_lr,
            maximal_learning_rate=max_lr,
            cycle_size=cycle_size,
            shift_peak=shift_peak,
            final_lr_scale=final_lr_scale
        )


# Nadam_beta_1
max_momentum = float(params['max_momentum'])
min_momentum = float(params['min_momentum'])

momentum_schedule = OneCycle(initial_learning_rate=max_momentum,
                             maximal_learning_rate=min_momentum,
                             cycle_size=cycle_size,
                             shift_peak=shift_peak,
                             final_lr_scale=1.0
                             )

# loss
huber_delta = float(params['huber_delta'])
loss_fxn = Huber(delta=huber_delta)

# optimizer:
# To implement OneCycle policy for LR and momentum, while also using weight decay
# legacy version of Nadam has to be used, as only this version has ._set_hyper method.
# However, legacy Nadam does not have weight decay, which is therefore extended using a tensorflow addon
NadamW = extend_with_decoupled_weight_decay(Nadam)
optimizer_fxn = NadamW(learning_rate=min_lr, # starting value of scheduler
                  beta_1=max_momentum, # starting value of scheduler
                  beta_2=Nadam_beta_2,
                  epsilon=Nadam_epsilon,
                  weight_decay=Nadam_weight_decay)

optimizer_fxn._set_hyper("learning_rate", lambda: LR_schedule(optimizer_fxn.iterations))
optimizer_fxn._set_hyper("beta_1", lambda: momentum_schedule(optimizer_fxn.iterations))


# custiom callback for reporting learning rate and momentum
if logging:
    class report_LR:
        def __init__(self, model, lrpath):
            self.model = model
            self.lrpath = lrpath
            self.file = open(lrpath, 'w')
            self.file.write("\t".join(["step", "lr", "beta_1\n"]))
            self.file.close()
        def on_batch_end(self, batch, logs):
            if batch % 10 == 0:
                self.file = open(self.lrpath, 'a')
                logs["lr"] = self.model.optimizer.lr
                logs["beta_1"] = self.model.optimizer.beta_1
                self.file.write("\t".join([str(batch), str(tf.get_static_value(logs["lr"])), str(tf.get_static_value(logs["beta_1"]))+"\n"]))
                self.file.close()
    rlr = report_LR(model, lrpath)
    rlr_cb = LambdaCallback(on_batch_end=lambda batch,
                            logs: rlr.on_batch_end(batch, logs))


model.compile(optimizer_fxn,
              loss = loss_fxn,
              metrics=[R2Score()])

if model_utils.to_bool(config['use_one_cycle_optimization']):
    if logging:
        cb_list.append(rlr_cb)
    num_epochs=cycle_epochs
else:
    cb_list.append(check_cb)
    cb_list.append(redLR_cb)
    cb_list.append(earlystop_cb)
    num_epochs=500 # maximum number of epochs

# train model
batch=int(float(params['batch_size']))

if full or (valid_fold is None):
    if full:
        print(f"Training model on the entire data set")
    else:
        print(f"Holding out fold {str(test_fold)} as a test set and training on the rest")
    print(f"Organism: {organism}")
    model.fit([train[i] for i in input_names],
              train["output"],
              batch_size=batch,
              epochs=num_epochs,
              validation_data=None,
              callbacks=cb_list,
              verbose=1) # display progress bars
else:
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


if model_utils.to_bool(config['use_one_cycle_optimization']):
    model.save_weights(os.path.join(outdir,model_outfile))
    #else saved via checkpoint cb
