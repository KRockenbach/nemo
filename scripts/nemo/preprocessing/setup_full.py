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
title: setup_full.py
description: sets up full training set for training nemo100, also fits and applies scalers and performs log transform
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      python nemo/preprocessing/setup_full.py <merged_data> <output_directory>
=========================================================================================================
'''

import sys, os
import pathlib
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn_pandas import DataFrameMapper
from pickle import dump
import math
import random

random.seed(1234)

##################################################
# AUXILLARY FUNCTION DEFINITIONS
##################################################

##################################################
# get entire log-transformed data table (so that it can remain in memory)
def get_data(data_file):
    '''
    loads data merged data and log-transforms appropriate columns
    '''
    # get entire data set to be held in memory for further processing
    print("+++++++++++++++++++++++++++++++++++++++++++++++++")
    table = pd.read_table(data_file, index_col=0) # read data, gene ID used as index.
    #index later used to select based on graphpart cluster

    table = table.loc[:,["MEDIAN_EXPRESSION","PROMOTER","TERMINATOR"]]
    table.iloc[:,0] = np.log10(table.iloc[:,0]+0.1)
    assert (not table.isnull().any().any()) # assert that table is free of NaN values
    return table


# get scaler and mapper for specific fold combination
def get_scaler(train_table, out_dir):
    num_out = 1
    # expression data for current fold
    exp = train_table.iloc[:,0:num_out]

    # sclaer to be saved for later inverse_transformation of prediction results
    # rest of numeric data is sclaed separately using mapper

    scaler = preprocessing.StandardScaler()

    # only train data is used to fit expression scaler
    scaler.fit(exp.values)
    # save scaler
    scaler_dir = os.path.join(out_dir, 'scalers')
    os.makedirs(scaler_dir, exist_ok=True) # equivalent to mkdir -p
    dump(scaler, open(os.path.join(scaler_dir, f'scaler_full.pkl'), 'wb'))
    return

#######################################


#######################################
# MAIN FUNCTION
#######################################

# import arguments
data_file = sys.argv[1]
out_dir = sys.argv[2]

# create output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


table = get_data(data_file)
# fit scaler and mapper
print('Total Number of training samples: {}'.format(table.shape[0]))
print('Fitting scaler')
get_scaler(table, out_dir)
print("Saving table")
table.to_feather(os.path.join(out_dir, f"full.feather"))

