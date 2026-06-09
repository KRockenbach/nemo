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
title: eval.py
description: makes predictions using trained model on corresponding test set
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-02-05
version: 1.0.2
usage: python -m nemo.eval.eval <directory containing data folds> <directory containing the model weights>
notes: run within nemo environment
=========================================================================================================
'''

import sys, h5py, os
import numpy as np
import pandas as pd
from random import choice
import csv
import tensorflow as tf
from pickle import load
from scipy import stats
from tensorflow.keras.models import load_model, Model
from tensorflow.keras.layers import Concatenate, Dense, Input
from sklearn.preprocessing import StandardScaler
from ..utils.model_utils import *
from ..models.nemo import build_nemo

def randomize_sequence(seq, segment="upstream", first_downstream_pos=5000):
    first_downstream_pos=int(first_downstream_pos)
    rand_seq = np.zeros_like(seq)
    for j in range(0,first_downstream_pos): # seq.shape[0]
        for i in range(seq.shape[1]):
            if seq[j,i,:].sum() == 1 and segment=="upstream" and j < first_downstream_pos:
                base = choice(range(4))
                rand_seq[j,i,base] = 1
            elif seq[j,i,:].sum() == 1 and segment=="downstream" and j >= first_downstream_pos:
                base = choice(range(4))
                rand_seq[j,i,base] = 1
            else:
                rand_seq[j,i,:] = seq[j,i,:]
    assert(rand_seq.sum() > 0)
    assert(rand_seq.sum() <= (seq.shape[0]*seq.shape[1]))
    return rand_seq


def inverse_transform(z, scaler):
        n_outputs=1
        z=z.reshape(-1,n_outputs)
        return (z*scaler.scale_[out_idx])+scaler.mean_[out_idx]


gpu_devices = tf.config.experimental.list_physical_devices('GPU')
for device in gpu_devices:
    tf.config.experimental.set_memory_growth(device, True)


folddir = "../data/Bnapus/masked_graphpart_fold_data"
weightdir = "../model_weights/nemo/Bnapus/masked_graphpart"
modelname = "nemo"
rootdir = ".."
modeldir = os.path.join(rootdir, "model_configs", modelname)
resultdir = modeldir.replace("model_configs", "results")
datadir = os.path.join(*folddir.split("/")[0:(len(folddir.split("/"))-1)])
test_descriptor = "masked_graphpart"
masking = "masked"
partitioning = "graphpart"
train_descriptor = masking + "_" + partitioning
organism = "Bnapus"
model_org = weightdir.split("/")[-2]
org_initial = model_org[0].upper()

outdir = os.path.join(resultdir, organism, train_descriptor, f"randomization")
os.makedirs(outdir, exist_ok=True)


# tsv file containing hyperparameters (only needed for batch size)
hp_path = os.path.join(modeldir, "hyperparams.tsv")


# read config
conf_path = os.path.join(modeldir, "config.tsv")
config = dict_from_tsv(conf_path)


# read hyperparams
params = dict_from_tsv(hp_path)
batch=int(float(params['batch_size']))
input_names = ["promoter", "terminator"]
outP = int(params['outside_P'])
inP = int(params['inside_P'])
outT = int(params['outside_T'])
inT = int(params['inside_T'])

predict_max = False
out_name = "median" # default
out_idx = 2

randomize=True

###################
def evaluate(test_fold, valid_fold, model_file, scaler, randomize, N):
    # load test data
    test = get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=folddir, set="test",
                  test_fold=test_fold, valid_fold=valid_fold)

    test["output"] = test["output"][:,out_idx]
    print(test["output"].shape)
    for key in test:
        exec(f'{key} = test["{key}"]') # output and ID are defined here

    model = build_nemo()

    model.load_weights(model_file)

    print('Loaded weights from:', model_file)

    # get list of input names
    input_names = list(test.keys())
    # get rid of output and ID
    input_names = input_names[:(len(input_names)-2)]
    print(input_names)

    # get list of input data
    inputs = []
    if randomize:
        rand_prom_up_inputs = []
        rand_term_up_inputs = []
        rand_prom_down_inputs = []
        rand_term_down_inputs = []
    for name in input_names:
        exec("inputs.append(" + name + ")")
        if name == "promoter" and randomize:
            rand_prom_up_inputs.append(randomize_sequence(test["promoter"], segment="upstream", first_downstream_pos=5000))
            rand_prom_down_inputs.append(randomize_sequence(test["promoter"], segment="downstream", first_downstream_pos=5000))
        elif randomize:
            exec("rand_prom_up_inputs.append(" + name + ")")
            exec("rand_prom_down_inputs.append(" + name + ")")
        if name == "terminator" and randomize:
            rand_term_up_inputs.append(randomize_sequence(test["terminator"], segment="upstream", first_downstream_pos=1200))
            rand_term_down_inputs.append(randomize_sequence(test["terminator"], segment="downstream", first_downstream_pos=1200))
        elif randomize:
            exec("rand_term_up_inputs.append(" + name + ")")
            exec("rand_term_down_inputs.append(" + name + ")")

    # get regression predictions
    #preds = model.predict(inputs, batch_size=batch)
    if randomize:
        rand_prom_up_preds = model.predict(rand_prom_up_inputs, batch_size=batch)
        rand_term_up_preds = model.predict(rand_term_up_inputs, batch_size=batch)
        rand_prom_down_preds = model.predict(rand_prom_down_inputs, batch_size=batch)
        rand_term_down_preds = model.predict(rand_term_down_inputs, batch_size=batch)

    ## perform inverse transformation
    x = inverse_transform(test["output"], scaler)
    #y = inverse_transform(preds, scaler)
    #x = scaler.inverse_transform(test["output"].reshape(-1,n_outputs))  #scaler expects 2D-array
    #y = scaler.inverse_transform(preds.reshape(-1,n_outputs))
    if randomize:
        rand_prom_up_y = inverse_transform(rand_prom_up_preds, scaler)
        rand_term_up_y = inverse_transform(rand_term_up_preds, scaler)
        rand_prom_down_y = inverse_transform(rand_prom_down_preds, scaler)
        rand_term_down_y = inverse_transform(rand_term_down_preds, scaler)
        #rand_prom_up_y = scaler.inverse_transform(rand_prom_up_preds.reshape(-1,n_outputs))
        #rand_term_up_y = scaler.inverse_transform(rand_term_up_preds.reshape(-1,n_outputs))
        #rand_prom_down_y = scaler.inverse_transform(rand_prom_down_preds.reshape(-1,n_outputs))
        #rand_term_down_y = scaler.inverse_transform(rand_term_down_preds.reshape(-1,n_outputs))


    gene_names = translate_IDs(test["ID"], datadir)
    #mat = np.column_stack((gene_names, x))

    colnames = ['Gene', 'Actual_Median', 'Predicted_Median']
    #if predict_max:
    #    colnames = ['Gene', 'Max_Expression']

    #df = pd.DataFrame(mat, columns=colnames)
    #f_out = os.path.join(outdir, f'actual.t_{test_fold}_{out_name}.txt')
    # actual expression only needs to be saved once per test fold
    #df.to_csv(f_out, index=False, header=True, sep='\t')
    #print(f"saved actual values to {f_out}")


    #mat = np.column_stack((gene_names, y))
    #df = pd.DataFrame(mat, columns=colnames)
    #f_out = os.path.join(outdir, f'predictions.t_{test_fold}.txt')
    #if valid_fold is not None:
    #    f_out = f_out.replace(".txt", f"_v_{valid_fold}.txt")
    #if N is not None:
    #    f_out = f_out.replace(".txt", f"_n_{N}_{out_name}.txt")
    # only one training rep per fold configuration, no further selection necessary
    #df.to_csv(f_out, index=False, header=True, sep='\t')
    #print(N)
    #print(valid_fold)
    #print(f"saved predicted values to {f_out}")

    if randomize: # N is not None
        mat = np.column_stack((gene_names, x, rand_prom_up_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_prom_up_predictions.t_{str(test_fold)}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, x, rand_prom_down_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_prom_down_predictions.t_{str(test_fold)}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, x, rand_term_up_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_term_up_predictions.t_{str(test_fold)}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, x, rand_term_down_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_term_down_predictions.t_{str(test_fold)}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")


###################

for test_fold in range(10):
    valid_fold = None
    model_file = os.path.join(weightdir, f"{modelname}_t_{str(test_fold)}_median.h5")
    scaler = load(open(os.path.join(folddir, 'scalers', f'scaler_{test_fold}.pkl'), 'rb'))
    randomize = True
    evaluate(test_fold, valid_fold, model_file, scaler, randomize, None)
