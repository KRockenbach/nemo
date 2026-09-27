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
date: 2026-07-29
version: 1.0.1
usage: python nemo/predict.py <directory containing data> <model file; defaults to Bnapus full model>
notes: run within nemo or nemo_cpu environment
=========================================================================================================
'''

import sys, os
import numpy as np
import pandas as pd
#from random import choice
#import csv
#import tensorflow as tf
from pickle import load
#from scipy import stats
from tensorflow.keras.models import load_model
from sklearn.preprocessing import StandardScaler
from models.nemo import build_nemo
#from utils.model_utils import *

###############################################################

def one_hot(seq):
    """ Takes in set of fasta sequences and converts them into
    3D array of one-hot encodings"""
    seq_len = len(seq.iloc[0])
    seqindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    seq_vec = np.zeros((seq.count(),seq_len,4), dtype='bool')
    for i in range(seq.count()):
        thisseq = seq.iloc[i]
        for j in range(seq_len):
            try:
                seq_vec[i,j,seqindex[thisseq[j]]] = 1
            except:
                # unknown nucleotides are indicated by 'N'
                # corresponding column is left as is, containing only zeros
                pass
    return seq_vec

def get_set4pred(datadir):
    data_dict = {}
    path = os.path.join(datadir, 'full.feather')
    table=pd.read_feather(path)
    data_dict["promoter"] = one_hot(table.loc[:,'PROMOTER'].str.slice(start=0,stop=6200))
    data_dict["terminator"] = one_hot(table.loc[:,'TERMINATOR'].str.slice(start=0,stop=6200))
    data_dict["ID"] = pd.Series(table.index).to_numpy(dtype=int)
    return(data_dict)

def inverse_transform(z, scaler):
        n_outputs=1
        z=z.reshape(-1,n_outputs)
        return (10**((z*scaler.scale_)+scaler.mean_) - 0.1)

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

##########################################################

def predict(datadir, model_file, scaler):
    modelname = "nemo"
    outdir = datadir
    batch=120

    # load test data
    test = get_set4pred(datadir=datadir)
    inputs = [test["promoter"], test["terminator"]]
    ID = test["ID"]

    model = build_nemo()
    model.load_weights(model_file)

    # get regression predictions
    preds = model.predict(inputs, batch_size=batch)

    gene_names = translate_IDs(test["ID"], datadir)

    y = inverse_transform(preds, scaler)
    mat = np.column_stack((gene_names, y))
    colnames = ["Gene", "Prediction"]
    df = pd.DataFrame(mat, columns=colnames)
    f_out = os.path.join(outdir, 'predictions.txt')
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(f"saved predicted values to {f_out}")

###########################################################

if __name__ == "__main__":
    datadir = sys.argv[1]
    model_file = sys.argv[2]
    scaler = load(open(sys.argv[3], 'rb'))
    predict(datadir, model_file, scaler)
