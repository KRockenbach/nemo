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
title: setup_folds.py
description: sets up training folds in feather format, also fits and applies scalers and performs log transform
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-08
version: 2.0.1
usage:
      python nemo/preprocessing/setup_folds.py <merged_data> <output_directory> <TSV with gene names and gene indices> <graphpart_output> <masking> <species>
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

num_outputs = 5

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

    table.dropna(axis=0, inplace=True) # remove any rows containing NaN values
    # expression output (num_outputs) and lengths (+4) and exon density (+1) 5UTR, 3UTR, CDS GC content (+3) promoter, terminator (+2)
    if len(table.columns) == (num_outputs+8+2):
        table.iloc[:,0:(num_outputs+5)] = np.log10(table.iloc[:,0:(num_outputs+5)]+0.1) # log_10 transformation
    else:
        table.iloc[:,0:num_outputs] = np.log10(table.iloc[:,0:num_outputs]+0.1) # log_10 transformation
    assert (not table.isnull().any().any())
    return table


# get fold-specific data
def get_fold_data(table, id_key, graphpart, fold, halflife_data=False):
    '''
    selects data for each fold based on graphpart results
    '''
    # for Bnapus, homoeologs with 100% identity have already been removed in prep_Bnapus.sh
    ########################################
    # get graphpart sets
    ########################################

    id_key.columns = ["gene_name", "ID"]
    # associate IDs in merged data table with gene names in graphpart table
    graphpart["gene_name"]=graphpart["AC"].apply(lambda x: x.split('.')[0]) # ensure that names are gene-level, not transcript-level
    merged_df = pd.merge(graphpart, id_key, how="inner", on="gene_name")
    gp_fold_idx = merged_df.loc[merged_df.loc[:,"cluster"].astype("int") == fold, "ID"].to_list()
    gp_fold_idx = [i for i in gp_fold_idx if i in table.index.to_list()] # NaNs were removed from table, need to be removed from index list too

    # data will be shuffled when loading it for training.
    # create data fold
    fold_table = table.copy().loc[gp_fold_idx,:]

    assert fold_table.index.to_list() == gp_fold_idx
    # print number of rows for verification
    if not halflife_data:
        if num_outputs == 1:
            fold_table = fold_table.loc[:,["MEDIAN_EXPRESSION", "PROMOTER","TERMINATOR"]]
        if num_outputs == 3:
            fold_table = fold_table.loc[:,["MIN_EXPRESSION", "MEDIAN_EXPRESSION", "MAX_EXPRESSION", "PROMOTER","TERMINATOR"]]
        else:
            fold_table = fold_table.loc[:,["MIN_EXPRESSION", "Q1_EXPRESSION", "MEDIAN_EXPRESSION", "Q3_EXPRESSION", "MAX_EXPRESSION", "PROMOTER","TERMINATOR"]]
    return fold_table


# get scaler and mapper for specific fold combination
def get_scaler(train_table, out_dir, test_fold, valid_fold):

    # expression data for current fold
    exp = train_table.iloc[:,0:num_outputs]

    # sclaer to be saved for later inverse_transformation of prediction results
    # rest of numeric data is sclaed separately using mapper

    scaler = preprocessing.StandardScaler()
    if len(train_table.columns) == (num_outputs+8+2):
        mapper = DataFrameMapper([(train_table.columns[0:num_outputs],
                                   None), # expression, scaled separately
                                  (train_table.columns[num_outputs:(num_outputs+8)], # 8 halflife features
                                   preprocessing.StandardScaler()), # scaled by mapper
                                  (train_table.columns[(num_outputs+8):],
                                   None)]) # promoter & terminator

    # only train data is used to fit expression scaler
    scaler.fit(exp.values)
    # save scaler
    scaler_dir = os.path.join(out_dir, 'scalers')
    os.makedirs(scaler_dir, exist_ok=True) # equivalent to mkdir -p
    if valid_fold is None:
        dump(scaler, open(os.path.join(scaler_dir, f'scaler_{test_fold}.pkl'), 'wb'))
        if len(train_table.columns) == 11:
            # fit on all numeric columns of train table
            mapper.fit(train_table)
            # save mapper
            dump(mapper, open(os.path.join(scaler_dir, f'mapper_{test_fold}.pkl'), 'wb'))
    else:
        dump(scaler, open(os.path.join(scaler_dir, f'scaler_{test_fold}_{valid_fold}.pkl'), 'wb'))
        if len(train_table.columns) == 11:
            # fit on all numeric columns of train table
            mapper.fit(train_table)
            # save mapper
            dump(mapper, open(os.path.join(scaler_dir, f'mapper_{test_fold}_{valid_fold}.pkl'), 'wb'))
    return

#######################################


#######################################
# MAIN FUNCTION
#######################################

# import arguments
data_file = sys.argv[1]
out_dir = sys.argv[2]
id_key = pd.read_table(sys.argv[3], index_col=False, header=None) # table containing gene names and corresponding IDs
gp = sys.argv[4]
graphpart = pd.read_table(gp, index_col=False, header=0, sep=',')

masking = sys.argv[5]
organism = sys.argv[6]


num_folds = 10


# create output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)


# load entire log-transformed data set into memory
# so it doesn't have to be loaded 90 times
table = get_data(data_file)

if masking == "masked" and organism == "Bnapus" and gp.split("/")[-1] == "graphpart_result.csv":
    Bn_gp = gp.replace(".csv", ".Bn.csv")
    Bn_graphpart = pd.read_table(Bn_gp, index_col=False, header=0, sep=',')
    Bn_out_dir=out_dir.replace("graphpart", "graphpart_Bn")
    if not os.path.exists(Bn_out_dir):
        os.makedirs(Bn_out_dir)
    test_fold=0
    for valid_fold in [1, None]:
        train_folds = [f for f in range(num_folds) if f not in [test_fold, valid_fold]]
        fold_df_list = []
        print(f"Retreiving train folds for test fold {test_fold} and valid fold {valid_fold}")
        for fidx, fold in enumerate(train_folds): #get train table for given fold combination
            fold_df_list.append(get_fold_data(table, id_key, Bn_graphpart, fold, halflife_data=True))
        print("Concatenating training folds together")
        train_table = pd.concat(fold_df_list)
        # fit and save scalers and mappers for current train table
        # print number of rows for verification
        print('Total Number of training samples: {}'.format(train_table.shape[0]))
        print('Fitting scaler for current fold configuration')
        get_scaler(train_table, (Bn_out_dir), test_fold, valid_fold)
        # free memory back up and foce creation of new table for next fold combination
        del train_table
        print("Train table deleted from memory\n")
    print("Saving folds")
    for fold in range(num_folds):
        fold_table = get_fold_data(table, id_key, Bn_graphpart, fold, halflife_data=True)
        fold_count = fold_table.shape[0]
        print(f'Number of samples in fold {fold}: {fold_count}')
        fold_table.to_feather(os.path.join(Bn_out_dir, f"fold_{fold}.feather"))


# fit scalers and mappers for all 10 possible combinations of test folds without valid folds (needed to train nemo with one-cycle learning rate policy)
for test_fold in range(num_folds):
    train_folds = [f for f in range(num_folds) if f != test_fold]
    fold_df_list = []
    print(f"Retreiving train folds for test fold {test_fold}")
    for fidx, fold in enumerate(train_folds): #get train table for given fold combination
        fold_df_list.append(get_fold_data(table, id_key, graphpart, fold))
    print("Concatenating training folds together")
    train_table = pd.concat(fold_df_list)
    # fit and save scalers and mappers for current train table
    # print number of rows for verification
    print('Total Number of training samples: {}'.format(train_table.shape[0]))
    print('Fitting scaler for current fold configuration')
    get_scaler(train_table, out_dir, test_fold, None)
    # free memory back up and foce creation of new table for next fold combination
    del train_table
    print("Train table deleted from memory\n")

print("Saving folds")
for fold in range(num_folds):
    fold_table = get_fold_data(table, id_key, graphpart, fold)
    fold_count = fold_table.shape[0]
    print(f'Number of samples in fold {fold}: {fold_count}')
    fold_table.to_feather(os.path.join(out_dir, f"fold_{fold}.feather"))


