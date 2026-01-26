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
title: model_utils.py
description: helper functions for data loading and model building
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-20
version: 1.0.1
usage:
      (1) run a given script (e.g. script.py) as: python -m script
      (2) within script, import functions using: from ..utils.model_utils import *
=========================================================================================================
'''

import h5py, os, sys, logging
logging.disable(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import random
import numpy as np
import pandas as pd
from pickle import load
import csv
from math import ceil, floor
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, AveragePooling1D, Flatten, Concatenate, Dense, BatchNormalization, Dropout, PReLU, AlphaDropout, ELU, Bidirectional, LSTM, Add, Multiply, LeakyReLU
from tensorflow.keras.activations import exponential, relu, selu, gelu, swish, softplus, mish
from tensorflow.keras.initializers import Constant
from tensorflow_addons.layers import GELU

tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)

num_out = 5

######################
def to_bool(string):
    if isinstance(string, bool): # check if input is already boolean
        return string
    elif string.upper() == "TRUE":
        return True
    elif string.upper() == "FALSE":
        return False
    else:
        raise ValueError

######################
def dict_from_tsv(tsv_path, verbose=False):
    reader = csv.reader(open(tsv_path, 'r'), delimiter='\t')
    dict = {}
    if verbose:
        print("Reading tsv:")
    for key, value in reader:
        dict[key] = value
        if verbose:
            print(key + ": " + value)
    if verbose:
        print("\n\n")
    return dict


######################
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


#####################
def translate_IDs(ID, datadir, verbose=False):
    # ID is a numpy array with dimensions (#samples, None)
    ID_series = pd.Series(ID[:,None].flatten())
    if verbose:
        print(f"ID shape: {str(ID_series.shape)}")
    ID_series.name = "ID"
    key_file = os.path.join(datadir, "derived_data", "gene.id.key")
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

######################
def load_fold(datadir, fold, scaler, mapper=None):
    # loads and scales data of a given fold, using appropriate scaler
    path = os.path.join(datadir, f'fold_{fold}.feather')
    fold_table=pd.read_feather(path)

    if mapper is not None:
        # scale half-life data (expression data remains unscaled for now)
        scaled_fold = mapper.transform(fold_table.copy())
        scaled_fold_table = pd.DataFrame(scaled_fold, index=fold_table.index,
                                         columns=fold_table.columns)

    else:
        scaled_fold_table = fold_table
    # scale expression data (done separately to ensure consistency with inverse transform)
    # unscaled expression data for current fold
    exp = scaled_fold_table.iloc[:,0:num_out]
    # scale expression data
    scaled_fold_table.iloc[:,0:num_out] = scaler.transform(exp.values)

    return scaled_fold_table


######################
def load_scaler(datadir, test_fold=0, valid_fold=1, verbose=False):
    if test_fold is None:
        scaler = load(open(os.path.join(datadir, 'scalers', f'scaler_full.pkl'), 'rb'))
        if verbose:
            print(f"Loaded scaler for the full set")
    elif valid_fold is None:
        scaler = load(open(os.path.join(datadir, 'scalers', f'scaler_{test_fold}.pkl'), 'rb'))
        if verbose:
            print(f"Loaded scaler for test fold {test_fold}")
    else:
        scaler = load(open(os.path.join(datadir, 'scalers', f'scaler_{test_fold}_{valid_fold}.pkl'), 'rb'))
        if verbose:
            print(f"Loaded scaler for test fold {test_fold} and valid fold {valid_fold}")
    return scaler


def load_mapper(datadir, test_fold=0, valid_fold=1, verbose=False):
    if test_fold is None:
        mapper = load(open(os.path.join(datadir, 'scalers', f'mapper_full.pkl'), 'rb'))
        if verbose:
            print(f"Loaded mapper for the full set")
    elif valid_fold is None:
        mapper = load(open(os.path.join(datadir, 'scalers', f'mapper_{test_fold}.pkl'), 'rb'))
        if verbose:
            print(f"Loaded mapper for test fold {test_fold}")
    else:
        mapper = load(open(os.path.join(datadir, 'scalers', f'mapper_{test_fold}_{valid_fold}.pkl'), 'rb'))
        if verbose:
            print(f"Loaded mapper for test fold {test_fold} and valid fold {valid_fold}")
    return mapper


######################
def retreive_column(scaled_fold_table, retreive="promoter", outP=15000, inP=5000, outT=15000, inT=5000, verbose=False):
    # retreive a specific column from the dataset
    # selection of columns depends on config and task
    if retreive == "promoter":
        full_len = scaled_fold_table.loc[:,'PROMOTER'].str.len().iloc[0]
        if verbose:
            print(f'Full {retreive}-proximal sequence length: {str(full_len)}')
            print(f'Number of examples: ' + str(scaled_fold_table.loc[:,'PROMOTER'].shape[0]))
        TSS=15000
        start=TSS-outP
        stop=TSS+inP
        if verbose:
            print("Subsampling and one-hot encoding promoter-proximal sequences...")
        promoters = one_hot(scaled_fold_table.loc[:,'PROMOTER'].str.slice(start=start, stop=stop))
        if verbose:
            print(f'Subsampled to [-{outP},{inP}] relative to TSS')
            print("New length = " + str(promoters.shape[1]))
        return promoters
    elif retreive == "terminator":
        full_len = scaled_fold_table.loc[:,'TERMINATOR'].str.len().iloc[0]
        if verbose:
            print(f'Full {retreive}-proximal sequence length: 20000')
            print(f'Number of examples: ' + str(scaled_fold_table.loc[:,'TERMINATOR'].shape[0]))
        TTS=5000
        start=TTS-inT
        stop=TTS+outT
        if verbose:
            print("Subsampling and one-hot encoding terminator-proximal sequences...")
        terminators = one_hot(scaled_fold_table.loc[:,'TERMINATOR'].str.slice(start=start, stop=stop))
        if verbose:
            print(f'Subsampled to: [-{inT},{outT}] relative to TTS')
            print("New length = " + str(terminators.shape[1]))
        return terminators
    elif retreive == "halflife":
        halflifedata = scaled_fold_table.iloc[:,num_out:(num_out+8)].to_numpy(dtype=float)
        return halflifedata
    elif retreive == "expression":
        expression = scaled_fold_table.iloc[:,0:num_out].to_numpy(dtype=float)
        return expression
    elif retreive == "ID":
        geneName = pd.Series(scaled_fold_table.index).to_numpy(dtype=int)
        return geneName



######################
def get_set(config, outP=15000, inP=5000, outT=15000, inT=5000, datadir="../data/Bnapus/masked_graphpart_fold_data", set="training", test_fold=0, valid_fold=1, verbose=False):
    # config file needed to determine which inputs to use
    input_names = []

    use_promoter = to_bool(config["use_promoter"])
    if use_promoter:
        input_names.append("promoter")
    use_terminator = to_bool(config["use_terminator"])
    if use_terminator:
        input_names.append("terminator")
    use_halflife = to_bool(config["use_halflife"])
    if use_halflife:
        input_names.append("halflife")

    if set == "full":
        folds = list(range(10))
        test_fold = None
        valid_fold= None
    elif set == "training" or set == "train":
        folds = [f for f in range(10) if f not in [test_fold, valid_fold]]
    elif set == "validation" or set == "valid":
        if valid_fold is None:
            return None
        else:
            folds = [valid_fold]
    elif set == "testing" or set == "test":
        folds = [test_fold]

    data_dict = {}

    scaler = load_scaler(datadir=datadir,
                         valid_fold=valid_fold,
                         test_fold=test_fold)

    if use_halflife:
        mapper = load_mapper(datadir=datadir,
                             valid_fold=valid_fold,
                             test_fold=test_fold)
    else:
        mapper = None

    first = True
    for fold in folds:
        fold_data = load_fold(datadir=datadir, fold=fold,
                              scaler=scaler, mapper=mapper)
        if verbose:
            print(f"\nLoading fold {fold}")

        for input in input_names:
            # load respective input columns
            input_data = retreive_column(fold_data, retreive=input, outP=outP, inP=inP, outT=outT, inT=inT, verbose=verbose)
            if first:
                data_dict[input] = input_data
            else:
                data_dict[input] = np.concatenate((data_dict[input], input_data), axis=0)

        output_data = retreive_column(fold_data, retreive="expression", verbose=verbose)
        if first:
            data_dict["output"] = output_data
        else:
            data_dict["output"] = np.concatenate((data_dict["output"], output_data), axis=0)

        ID_data = retreive_column(fold_data, retreive="ID", verbose=verbose)
        if first:
            data_dict["ID"] = ID_data
        else:
            data_dict["ID"] = np.concatenate((data_dict["ID"], ID_data), axis=0)
        # switch flag for consecutive folds
        first=False
        if set=="full" and verbose :
            print(f"Data fully loaded\n")
        elif verbose :
            print(f"Data for fold {fold} fully loaded\n")

    # pseudo-randomly shuffle data with seed
    rng = random.Random(42)
    key_list = [key for key in data_dict]
    shuff_idx = list(range(len(data_dict[key_list[0]])))

    # shuffle indices inplace
    rng.shuffle(shuff_idx)

    shuff_data_dict = {}
    if verbose:
        print(f"Randomly shuffling {set} set")
    for key in key_list:
        sample_axis=0
        shuff_data_dict[key] = np.take(data_dict[key], shuff_idx, sample_axis)

    n_folds = len(folds)
    if verbose:
        print(f"{set} set fully loaded... ({n_folds} fold(s))\n")
    return(shuff_data_dict)



######################
def get_set4pred(datadir="../validation/ZS11/data"):

    data_dict = {}
    path = os.path.join(datadir, 'full.feather')
    table=pd.read_feather(path)

    data_dict["promoter"] = one_hot(table.loc[:,'PROMOTER'].str.slice(start=0,stop=6200))
    data_dict["terminator"] = one_hot(table.loc[:,'TERMINATOR'].str.slice(start=0,stop=6200))
    data_dict["ID"] = pd.Series(table.index).to_numpy(dtype=int)

    return(data_dict)


######################
def build_model(params, train, input_names=["promoter", "terminator"]):

    # Define Hyperparameters
    if "promoter" in input_names:
        ncl_P=int(float(params['ncl_P'])) # number of conv layers in promoter branch
        nlstm_P=int(float(params['nlstm_P']))
    if "terminator" in input_names:
        ncl_T=int(float(params['ncl_T'])) # number of conv layers in terminator branch
        nlstm_T=int(float(params['nlstm_T']))
    try:
        use_cb_lstm = to_bool(params['use_cb_lstm'])
    except KeyError:
        use_cb_lstm = False
    ndl=int(float(params['ndl'])) # number of regression dense layers
    use_Bnorm = to_bool(params['use_Bnorm'])
    if use_Bnorm: # quite resource intensive
        print("Using batch normalization.")
    use_drop = to_bool(params['use_drop'])
    if use_drop:
        use_conv_drop = to_bool(params['use_conv_drop'])
        if use_conv_drop:
            print("Using dropout after conv layer")
        else:
            print("Not using dropout after conv layer")
        try:
            use_lstm_drop = to_bool(params['use_lstm_drop'])
            if use_lstm_drop:
                print("Using dropout after LSTM")
            else:
                print("Not using dropout after LSTM")
        except KeyError:
            use_lstm_drop=False
    use_exp = to_bool(params['use_exp'])

    #define activation function and activation-specific initialization and dropout
    act_fxn = params["activation"]
    kernel_ini = str(params[f'kernel_ini_{act_fxn}']) # kernel initiation

    if act_fxn == "selu":
        # initialization must be LeCun_init
        # dropout must be alpha dropout
        def activation(x):
            x = selu(x)
            return(x)
        def drop(x, rate):
            x = AlphaDropout(rate=rate)(x)
            return(x)
    else:
        def drop(x, rate):
            x = Dropout(rate=rate)(x)
            return(x)
        if act_fxn == "leaky_relu" or act_fxn == "prelu":
            # initialization must be He_init
            relu_alpha = float(params['relu_alpha'])
            if act_fxn == "leaky_relu":
                def activation(x):
                    x = relu(x, alpha = relu_alpha)
                    return(x)
            else: # prelu
                # relu alpha serves as alpha initialization
                def activation(x):
                    x = PReLU(alpha_initializer=Constant(relu_alpha))(x)
                    return(x)
        elif act_fxn == "gelu":
            def activation(x):
                # in case of bad performance, it's possible to approximate
                x = gelu(x)
                return(x)
        elif act_fxn == "swish":
            def activation(x):
                x = swish(x)
                return(x)
        elif act_fxn == "mish":
            def activation(x):
                x = mish(x)
                return(x)
        elif act_fxn == "softplus":
            def activation(x):
                x = softplus(x)
                return(x)
        elif act_fxn == "elu":
            elu_alpha = float(params['elu_alpha'])
            def activation(x):
                x = ELU(alpha=elu_alpha)(x)
                return(x)
        else: #standard relu
            def activation(x):
                x = relu(x, alpha = 0.0)
                return(x)
    print("Activation: " + act_fxn)
    print("Initialization: " + kernel_ini)

    ###################################################################################
    # define model
    ###################################################################################
    inputs = {}
    branches = []
    for name in input_names:
        input = Input(shape=train[name].shape[1:], name=name)
        inputs[name]=input
        branch = name[0].upper() # P, T and H denote respective input branches
        exec(f'{branch} = input')
        # construction of convolutional layers
        if name in ["promoter", "terminator"]:
            ncl = eval(f'ncl_{branch}')
            print (f"Constructing {ncl} convolutional layers in {name} branch.")
            ########### construct convolutional layers
            for layer in range(ncl):
                pl_size = f'pool_size_{branch}_{layer}of{ncl}'
                pl_type = f'pool_type_{branch}_{layer}of{ncl}'
                flt = f'filter_{branch}_{layer}of{ncl}'
                krn = f'kernel_{branch}_{layer}of{ncl}'
                cv_frac = f'conv_frac_{branch}_{layer}of{ncl}'
                try:
                    cv_stride = max([1, ceil(float(params[cv_frac])*int(params[krn]))])
                except KeyError:
                    if layer == 0:
                        print("No 'conv_frac' entry found")
                    cv_stride = 1
                dilation_rate = 1 # used hardware doesn't permit dilated convolution
                # convolution
                exec(f'''{branch} = Conv1D(int(params[flt]),
                       int(params[krn]),
                       dilation_rate = dilation_rate,
                       strides = cv_stride,
                       padding = 'same',
                       kernel_initializer = kernel_ini)({branch})''')

                # batch normalization
                if use_Bnorm:
                    exec(f'{branch} = BatchNormalization(momentum=float(params["Bnorm_momentum"]))({branch})')
                # activation
                if use_exp and layer==0: # only used after first conv layer, if at all
                    exec(f'{branch} = exponential({branch})')
                    print("Using exponential activation after first layer!")
                else:
                    exec(f'{branch} = activation({branch})')
                # convolutional dropout
                if use_drop:
                    if use_conv_drop:
                        drop_rate = f'cv_drop_{branch}_{layer}of{ncl}'
                        exec(f'{branch} = drop({branch}, rate=float(params[drop_rate]))')
                # pooling
                if params[pl_type] == "avrg" and int(params[pl_size]) > 1:
                    pl_frac = f'pool_frac_{branch}_{layer}of{ncl}'
                    pl_stride = max([1, ceil(float(params[pl_frac])*int(params[pl_size]))])
                    exec(f'''{branch} = AveragePooling1D(int(params[pl_size]),
                                                 padding = 'same',
                                                 strides=pl_stride)({branch})''')
                #elif params[pl_type] == "attention":
                #    exec(f'''{branch} = SoftmaxPooling1D(pool_size=int(params[pl_size]), per_channel=True, w_init_scale=2.0)({branch})''')
                elif int(params[pl_size]) > 1:
                    pl_frac = f'pool_frac_{branch}_{layer}of{ncl}'
                    pl_stride = max([1, ceil(float(params[pl_frac])*int(params[pl_size]))])
                    exec(f'''{branch} = MaxPooling1D(int(params[pl_size]),
                                                 padding = 'same',
                                                 strides=pl_stride)({branch})''')
                else:
                    print(f"No pooling in branch {branch} layer {layer}")

            # lstm block
            nlstm = eval(f'nlstm_{branch}')
            if nlstm > 0:
                skip_type = params['skip_type']
                if skip_type != 'none':
                    print('using skip connection')
                    exec(f'skip_{branch} = {branch}')
                    if skip_type in ['add', 'multiply']:
                        exec(f'skip_{branch} = Concatenate(axis = 2)([skip_{branch}, skip_{branch}])') # one copy of the feature map for each direction of the bidirectional LSTM
                    exec(f'skip_{branch} = BatchNormalization(momentum=float(params["Bnorm_momentum"]))(skip_{branch})')
                for layer in range(nlstm):
                    exec(f'{branch} = Bidirectional(LSTM(int(params["lstm_units_{branch}_{layer}of{nlstm}"]), return_sequences=True))({branch})')
                    # batch normalization
                    if use_Bnorm:
                        exec(f'{branch} = BatchNormalization(momentum=float(params["Bnorm_momentum"]))({branch})')
                    # lstm dropout
                    if use_drop:
                        if use_lstm_drop:
                            exec(f'{branch} = Dropout(float(params["lstm_drop_{branch}_{layer}of{nlstm}"]))({branch})')
                if skip_type == 'add':
                    exec(f'{branch} = Add()([{branch}, skip_{branch}])')
                elif skip_type == 'multiply':
                    exec(f'{branch} = Multiply()([{branch}, skip_{branch}])')
                elif skip_type == 'concatenate':
                    exec(f'{branch} = Concatenate(axis=2)([{branch}, skip_{branch}])') # axis2 = units/filters
            if use_cb_lstm == False:
                # flatten if fed directly into dense layers
                exec(f'{branch} = Flatten()({branch})')
            exec(f'branches.append({branch})') # add convolutional branches

    if use_cb_lstm:
        print('using cross-branch LSTM')
        L = Concatenate(axis=1)(branches)
        L = Bidirectional(LSTM(int(params["cb_lstm_units"]), return_sequences=True))(L)
        # flatten lstm output, to be fed into dense layers
        D = Flatten()(L)
        if 'halflife' in input_names:
            D = Concatenate(axis=1)([D,H]) # concat halflife features
    else:
        # add halflife and concat
        if 'halflife' in input_names:
            branches.append(inputs['halflife'])
        D = Concatenate(axis=1)(branches)

    # define dense layers
    print(f"Constructing {ndl} dense layers.")
    for layer in range(ndl):
        dl_dim = f'dl_dim_{layer}of{ndl}'
        dl_drop = f'dl_drop_{layer}of{ndl}'
        D = Dense(int(params[dl_dim]), kernel_initializer=kernel_ini)(D)
        if use_Bnorm:
            D = BatchNormalization(momentum=float(params["Bnorm_momentum"]))(D)
        D = activation(D)
        if use_drop and float(params[dl_drop]) > 0.0:
            D = drop(D, float(params[dl_drop]))
    output = Dense(1, name='out')(D)

    model = Model(inputs=[inputs[name] for name in input_names],
                  outputs=output)

    print("Model has been fully constructed.")
    p_count = sum([layer.count_params() for layer in model.layers])
    print("Total number of model parameters: " + str(p_count))
    return model


######################################



def build_nemo_promoter(promoter, momentum=0.81669, pooling_fixed=False):

    ini = "glorot_normal"
    def get_stride(fraction, size):
        return max([1, ceil(fraction*size)])
    m = momentum

    P = Conv1D(288, 5, padding = 'same', kernel_initializer = ini)(promoter)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(106, 8, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)
    if pooling_fixed:
        P = AveragePooling1D(2*get_stride(0.7, 21), strides = get_stride(0.7, 21), padding="same")(P)
    else:
        P = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(P)

    P = Conv1D(227, 9, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(187, 10, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(64, 46, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)
    if pooling_fixed:
        P = AveragePooling1D(2*get_stride(0.7, 18), strides = get_stride(0.7, 18), padding="same")(P)
    else:
        P = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(P)

    P = Conv1D(152, 62, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=m)(P)
    P = GELU(approximate=False)(P)
    if pooling_fixed:
        P = AveragePooling1D(2*get_stride(0.7, 8), strides = get_stride(0.7, 8), padding="same")(P)
    else:
        P = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(P)

    P = Flatten()(P)
    return P

def build_nemo_terminator(terminator, momentum=0.81669, pooling_fixed=False):

    ini = "glorot_normal"
    def get_stride(fraction, size):
        return max([1, ceil(fraction*size)])
    m = momentum

    T = Conv1D(219, 4, padding = 'same', kernel_initializer = ini)(terminator)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(449, 4, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)
    if pooling_fixed:
        T = AveragePooling1D(2*get_stride(0.7, 21), strides = get_stride(0.7, 21), padding="same")(T)
    else:
        T = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(T)

    T = Conv1D(259, 16, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(154, 15, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(105, 20, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)
    if pooling_fixed:
        T = AveragePooling1D(2*get_stride(0.7, 18), strides = get_stride(0.7, 18), padding="same")(T)
    else:
        T = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(T)

    T = Conv1D(152, 22, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=m)(T)
    T = GELU(approximate=False)(T)
    if pooling_fixed:
        T = AveragePooling1D(2*get_stride(0.7, 8), strides = get_stride(0.7, 8), padding="same")(T)
    else:
        T = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(T)

    T = Flatten()(T)
    return T


def build_nemo():
    promoter = Input(shape=(6200, 4), name="promoter")
    terminator = Input(shape=(6200, 4), name="terminator")
    P = build_nemo_promoter(promoter)
    T = build_nemo_temrinator(terminator)
    D = Concatenate(axis=1)([P,T])
    D = Dense(750, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)

    D = Dense(3, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)
    D = Dropout(0.00124)(D)
    D = Dense(1)(D) # single-regression output
    return Model(inputs = [promoter, terminator], outputs = D)



def build_simple_encoder(input, momentum=0.81669):
    # pooling windows are fixed to double the stride
    # for "same" padding output size of pooling layer only depends on stride
    ini = "glorot_normal"
    pl_stride = 10

    m = momentum
    def conv_block(x, filter, kernel):
        x = Conv1D(filter, kernel, padding = 'same', kernel_initializer = ini)(x)
        x = BatchNormalization(momentum=m)(x)
        x = GELU(approximate=False)(x)

        # double the filters and double the kernel size
        x = Conv1D(2*filter, 2*kernel, padding = 'same', kernel_initializer = ini)(x)
        x = BatchNormalization(momentum=m)(x)
        x = GELU(approximate=False)(x)
        x = AveragePooling1D(2*pl_stride, strides = pl_stride, padding='same')(x)
        return x

    X = conv_block(input, 50, 5)
    X = conv_block(X, 100, 10)
    X = conv_block(X, 200, 31) # last conv kernels cover entire input


    X = Flatten()(X)
    return X

def build_nemo2_dense(P, T, dim_0, dim_1, dropout, momentum=0.81669):
    m = momentum
    ini = "glorot_normal"
    x = Concatenate(axis=1)([P, T])

    x = Dense(dim_0, kernel_initializer=ini)(x)
    x = BatchNormalization(momentum=m)(x)
    x = GELU(approximate=False)(x)

    x = Dense(dim_1, kernel_initializer=ini)(x)
    x = BatchNormalization(momentum=m)(x)
    x = GELU(approximate=False)(x)

    x = Dropout(dropout)(x)

    min_out = Dense(1, name="min_out")(x)
    Q1_out = Dense(1, name="Q1_out")(x)
    Q2_out = Dense(1, name="Q2_out")(x)
    Q3_out = Dense(1, name="Q3_out")(x)
    max_out = Dense(1, name="max_out")(x)
    outputs = [min_out, Q1_out, Q2_out, Q3_out, max_out]
    return outputs

def build_nemo2():
    promoter = Input(shape=(6200, 4), name="promoter")
    terminator = Input(shape=(6200, 4), name="terminator")
    P = build_nemo2_promoter_branch(...)
    T= build_nemo2_terminator_branch(...)
    D = Concatenate(axis=1)([P, T])
    D = Dense(750, kernel_initializer=ini)(encoder)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)

    D = Dense(3, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)
    D = Dropout(0.00124)(D)
    D = Dense(5)(D) # multi-regression output
    return Model(inputs = [promoter, terminator], outputs = D)
