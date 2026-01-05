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
date: 2026-01-05
version: 1.0.1
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
from tensorflow.keras.models import load_model
from sklearn.preprocessing import StandardScaler
from ..utils.model_utils import *


def randomize_sequence(seq, segment="upstream", first_downstream_pos="5000"):
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


gpu_devices = tf.config.experimental.list_physical_devices('GPU')
for device in gpu_devices:
    tf.config.experimental.set_memory_growth(device, True)


folddir = sys.argv[1]
weightdir = sys.argv[2]
modelname = weightdir.split("/")[-3]
rootdir = weightdir.split("/")[0]
modeldir = os.path.join(rootdir, "model_configs", modelname)
resultdir = modeldir.replace("model_configs", "results")
datadir = os.path.join(*folddir.split("/")[0:(len(folddir.split("/"))-1)])
test_descriptor = folddir.split("/")[-1].replace("_fold_data", "")
masking = test_descriptor.replace("_graphpart_Bn", "").replace("_graphpartDS", "").replace("_graphpart", "").replace("_random", "")
partitioning = test_descriptor.replace(f"{masking}_", "")
train_descriptor = weightdir.split("/")[-1]
organism = folddir.split("/")[-2]
model_org = weightdir.split("/")[-2]
org_initial = model_org[0].upper()

print(partitioning)
if partitioning == "graphpartDS":
    outdir = os.path.join(resultdir, organism, train_descriptor, f"{modelname}{org_initial}_DS_preds")
else:
    outdir = os.path.join(resultdir, organism, train_descriptor, f"{modelname}{org_initial}_preds")
os.makedirs(outdir, exist_ok=True)


# tsv file containing hyperparameters (only needed for batch size)
hp_path = os.path.join(modeldir, "hyperparams.tsv")


# read config
conf_path = os.path.join(modeldir, "config.tsv")
config = dict_from_tsv(conf_path)


# read hyperparams
params = dict_from_tsv(hp_path)
batch=int(float(params['batch_size']))
input_names = []
if to_bool(config["use_promoter"]):
    outP = int(params['outside_P'])
    inP = int(params['inside_P'])
    input_names.append("promoter")
else:
    outP = 0
    inP = 0
if to_bool(config["use_terminator"]):
    outT = int(params['outside_T'])
    inT = int(params['inside_T'])
    input_names.append("terminator")
else:
    outT = 0
    inT = 0
if to_bool(config["use_halflife"]):
    input_names.append("halflife")



###################
def evaluate(test_fold, valid_fold, model_file, scaler, randomize, N):
    # load test data
    test = get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=folddir, set="test",
                  test_fold=test_fold, valid_fold=valid_fold)
    for key in test:
        exec(f'{key} = test["{key}"]') # output and ID are defined here


    if modelname == "xpresso":
        model = build_xpresso()
    else:
        model = build_nemo()

    model.load_weights(model_file)

    print('Loaded weights from:', model_file)

    # get list of input names
    input_names = list(test.keys())
    # get rid of output and ID
    input_names = input_names[:(len(input_names)-2)]

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
            rand_prom_up_inputs.append(randomize_sequence(test["promoter"], segment="upstream", first_downstream_pos="5000"))
            rand_prom_down_inputs.append(randomize_sequence(test["promoter"], segment="downstream", first_downstream_pos="5000"))
        elif randomize:
            exec("rand_prom_up_inputs.append(" + name + ")")
            exec("rand_prom_down_inputs.append(" + name + ")")
        if name == "terminator" and randomize:
            rand_term_up_inputs.append(randomize_sequence(test["terminator"], segment="upstream", first_downstream_pos="1200"))
            rand_term_down_inputs.append(randomize_sequence(test["terminator"], segment="downstream", first_downstream_pos="1200"))
        elif randomize:
            exec("rand_term_up_inputs.append(" + name + ")")
            exec("rand_term_down_inputs.append(" + name + ")")

    # get regression predictions
    preds = model.predict(inputs, batch_size=batch)
    if randomize:
        rand_prom_up_preds = model.predict(rand_prom_up_inputs, batch_size=batch)
        rand_term_up_preds = model.predict(rand_term_up_inputs, batch_size=batch)
        rand_prom_down_preds = model.predict(rand_prom_down_inputs, batch_size=batch)
        rand_term_down_preds = model.predict(rand_term_down_inputs, batch_size=batch)

    # perform inverse transformation
    n_outputs = 3
    x = scaler.inverse_transform(test["output"].reshape(-1,n_outputs))  #scaler expects 2D-array
    y = scaler.inverse_transform(preds.reshape(-1,n_outputs))
    if randomize:
        rand_prom_up_y = scaler.inverse_transform(rand_prom_up_preds.reshape(-1,n_outputs))
        rand_term_up_y = scaler.inverse_transform(rand_term_up_preds.reshape(-1,n_outputs))
        rand_prom_down_y = scaler.inverse_transform(rand_prom_down_preds.reshape(-1,n_outputs))
        rand_term_down_y = scaler.inverse_transform(rand_term_down_preds.reshape(-1,n_outputs))


    gene_names = translate_IDs(test["ID"], datadir)
    mat = np.column_stack((gene_names, x))
    colnames = ['Gene', 'Median_Expression']
    df = pd.DataFrame(mat, columns=colnames)
    f_out = os.path.join(outdir, f'actual.t_{test_fold}.txt')
    # actual expression only needs to be saved once per test fold
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(f"saved actual values to {f_out}")


    mat = np.column_stack((gene_names, y))
    df = pd.DataFrame(mat, columns=colnames)
    f_out = os.path.join(outdir, f'predictions.t_{test_fold}.txt')
    if valid_fold is not None:
        f_out = f_out.replace(".txt", f"_v_{valid_fold}.txt")
    if N is not None:
        f_out = f_out.replace(".txt", f"_n_{N}.txt")
    # only one training rep per fold configuration, no further selection necessary
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(N)
    print(valid_fold)
    print(f"saved predicted values to {f_out}")

    if randomize: # N is not None
        mat = np.column_stack((gene_names, rand_prom_up_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_prom_up_predictions.t_{str(test_fold)}_v_{valid_fold}_n_{N}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, rand_prom_down_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_prom_down_predictions.t_{str(test_fold)}_v_{valid_fold}_n_{N}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, rand_term_up_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_term_up_predictions.t_{str(test_fold)}_v_{valid_fold}_n_{N}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")

        mat = np.column_stack((gene_names, rand_term_down_y))
        df = pd.DataFrame(mat, columns=colnames)
        f_out = os.path.join(outdir, f'rand_term_down_predictions.t_{str(test_fold)}_v_{valid_fold}_n_{N}.txt')
        # only one training rep per fold configuration, no further selection necessary
        df.to_csv(f_out, index=False, header=True, sep='\t')
        print(f"saved predicted values to {f_out}")


###################

for test_fold in range(10):
    if (partitioning == "graphpart_Bn" or modelname == "xpresso") and (test_fold != 0):
        break
    elif partitioning == "graphpart_Bn" or modelname == "xpresso":
        for N in range(10):
            valid_fold = 1
            model_file = os.path.join(weightdir, f"{modelname}_t_0_v_1_n_{N}.h5")
            scaler = load(open(os.path.join(folddir, 'scalers', 'scaler_0_1.pkl'), 'rb'))
            if modelname == "nemo":
                randomize = True
            else:
                randomize = False
            evaluate(test_fold, valid_fold, model_file, scaler, randomize, N)

            if modelname == "nemo":
                valid_fold = None
                model_file = os.path.join(weightdir, f"nemo_t_0_n_{N}.h5")
                scaler = load(open(os.path.join(folddir, 'scalers', 'scaler_0.pkl'), 'rb'))
                randomize = False
                evaluate(test_fold, valid_fold, model_file, scaler, randomize, N)
    else:
        valid_fold = None
        model_file = os.path.join(weightdir, f"{modelname}_t_{str(test_fold)}.h5")
        scaler = load(open(os.path.join(folddir, 'scalers', f'scaler_{test_fold}.pkl'), 'rb'))
        randomize = False
        evaluate(test_fold, valid_fold, model_file, scaler, randomize, None)





###########
# Get preds on full set

if partitioning == "graphpart" and masking == "masked" and modelname == "nemo":

    test = get_set(config, outP=outP, inP=inP, outT=outT, inT=inT, datadir=folddir, set="full",
                      test_fold=None, valid_fold=None)
    for key in test:
        exec(f'{key} = test["{key}"]') # output and ID are defined here


    #get best learned parameters
    model_file = os.path.join(weightdir, f"{modelname}_full.h5")
    model = build_nemo()
    model.load_weights(model_file)

    print('Loaded weights from:', model_file)

    # get list of input names
    input_names = list(test.keys())
    # get rid of output and ID
    input_names = input_names[:(len(input_names)-2)]
    # get list of input data
    inputs = []
    for name in input_names:
        exec("inputs.append(" + name + ")")
    # get regression predictions
    preds = model.predict(inputs, batch_size=batch)


    # get scaler for inverse transformation
    scaler = load(open(os.path.join(folddir, 'scalers', f'scaler_full.pkl'), 'rb'))

    # perform inverse transformation
    print(test["output"].shape)
    n_outputs = 3
    x = scaler.inverse_transform(test["output"].reshape(-1,n_outputs))  #scaler expects 2D-array
    y = scaler.inverse_transform(preds.reshape(-1,n_outputs))


    gene_names = translate_IDs(test["ID"], datadir)
    mat = np.column_stack((gene_names, x))
    colnames = ['Gene', 'Median_Expression']
    df = pd.DataFrame(mat, columns=colnames)
    f_out = os.path.join(outdir, f'actual.full.txt')
    # actual expression only needs to be saved once per test fold
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(f"saved actual values to {f_out}")


    mat = np.column_stack((gene_names, y))
    df = pd.DataFrame(mat, columns=colnames)
    f_out = os.path.join(outdir, f'predictions.full.txt')
    # only one training rep per fold configuration, no further selection necessary
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(f"saved predicted values to {f_out}")
