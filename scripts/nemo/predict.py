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
title: predict.py
description: makes predictions using trained model
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python nemo/predict.py <directory containing data> <model file; defaults to Bnapus full model>
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
from utils.model_utils import *


gpu_devices = tf.config.experimental.list_physical_devices('GPU')
for device in gpu_devices:
    tf.config.experimental.set_memory_growth(device, True)

root=".."

datadir = sys.argv[1]
if len(sys.argv) == 3:
    model_file = sys.argv[2]
else:
    model_file = f"{root}/model_weights/nemo/Bnapus/masked_graphpart/nemo_full.h5"

modelname = "nemo" #weightdir.split("/")[-3]
outdir = datadir #modeldir.replace("model_configs", "results")
batch=120


scaler = load(open(model_file.replace('nemo_full.h5', 'scaler_full.pkl'), 'rb'))

# load test data
test = get_set4pred(datadir=datadir)
inputs = [test["promoter"], test["terminator"]]
ID = test ["ID"]

model = build_nemo()
model.load_weights(model_file)
print('Loaded weights from:', model_file)

# get regression predictions
preds = model.predict(inputs, batch_size=batch)

def translate_IDs(ID, datadir, verbose=False):
    # ID is a numpy array with dimensions (#samples, None)
    ID_series = pd.Series(ID[:,None].flatten())
    if verbose:
        print(f"ID shape: {str(ID_series.shape)}")
    ID_series.name = "ID"
    key_file = os.path.join(datadir, "gene.id.key")
    key_df = pd.read_table(key_file, index_col=False, header=None)
    if verbose:
        print(f"key_df shape: {str(key_df.shape)}")
    key_df.columns = ["gene_name", "ID"]
    for ID in ID_series:
        if ID not in key_df.loc[:,"ID"]:
            print(ID)
    # do inner join, preserving the order of the IDs
    merged_df = pd.merge(ID_series, key_df, how="inner", on="ID")
    if verbose:
        print(f"merged_df shape: {str(merged_df.shape)}")
    gene_names = merged_df.loc[:,"gene_name"].to_list()
    return gene_names


gene_names = translate_IDs(test["ID"], datadir)



y = scaler.inverse_transform(preds.reshape(-1,1))
mat = np.column_stack((gene_names, y))
colnames = ["Gene", "Prediction"]
df = pd.DataFrame(mat, columns=colnames)
f_out = os.path.join(outdir, 'predictions.txt')
df.to_csv(f_out, index=False, header=True, sep='\t')
print(f"saved predicted values to {f_out}")

