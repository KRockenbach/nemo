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
title: attribution.py
description: calculates attributions for nemo90 models on masked (grahpart partitioned) sequences on respective test fold
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.attribution.attribution <species> <test_fold_number>
notes: run within nemo environment
=========================================================================================================
'''

import sys, h5py, os, logging
logging.disable(logging.WARNING)
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
import numpy as np
import pandas as pd
from tqdm.contrib import itertools # for progress bar
import tensorflow as tf
from tensorflow.keras import Model
from ..utils.model_utils import *
from ..utils.attrib_utils import *

tf.compat.v1.logging.set_verbosity(tf.compat.v1.logging.ERROR)


# random state for dinuc shuffle
r = np.random.RandomState(1234)

root = '..'

organism = sys.argv[1] # directory containing training data
datadir = os.path.join(root, 'data', f'{organism}')
model_name = "nemo"
modeldir = os.path.join(root, 'model_configs', model_name)
resultdir = modeldir.replace('model_configs', 'results')
data_descriptor = "masked_graphpart"
folddir = os.path.join(datadir, f'{data_descriptor}_fold_data')
conf_path = os.path.join(modeldir, 'config.tsv')
config = dict_from_tsv(conf_path)
outdir = os.path.join(resultdir, organism, 'masked_graphpart', 'attribs')
os.makedirs(outdir, exist_ok=True)
test_fold = int(sys.argv[2])
param_path = os.path.join(modeldir, 'hyperparams.tsv')
params = dict_from_tsv(param_path)

if to_bool(config['use_promoter']):
    outP = int(params['outside_P'])
    inP = int(params['inside_P'])
if to_bool(config['use_terminator']):
    outT = int(params['outside_T'])
    inT = int(params['inside_T'])


################
valid_fold=None
if to_bool(config['use_promoter']) and to_bool(config['use_terminator']):
    test = get_set(config, outP=outP, inP=inP, outT=outT, inT=inT,
                   datadir=folddir, set='test', test_fold=test_fold, valid_fold=valid_fold)
elif to_bool(config['use_promoter']) and not to_bool(config['use_terminator']):
    test = get_set(config, outP=outP, inP=inP,
                   datadir=folddir, set='test', test_fold=test_fold, valid_fold=valid_fold)
elif to_bool(config['use_terminator']) and not to_bool(config['use_promoter']):
    test = get_set(config, outT=outT, inT=inT,
                   datadir=folddir, set='test', test_fold=test_fold, valid_fold=valid_fold)

gene_names = translate_IDs(test['ID'], datadir=datadir)

inputs = []
input_types = []
input_names = []
if to_bool(config['use_promoter']):
    inputs.append(test['promoter'].astype('int'))
    input_types.append('seq')
    input_names.append('promoter')
if to_bool(config['use_terminator']):
    inputs.append(test['terminator'].astype('int'))
    input_types.append('seq')
    input_names.append('terminator')
if to_bool(config['use_halflife']):
    inputs.append(test['halflife'].astype('float'))
    input_types.append('num')
    input_names.append('halflife')
assert len(inputs) > 0

model_path = os.path.join(resultdir.replace("results", "model_weights"),
                          organism,
                          'masked_graphpart',
                          f'{model_name}_t_{str(test_fold)}.h5')
exec(f"model = build_{model_name}()")
model.load_weights(model_path)
attr = gradshap(model=model, x=inputs, input_type=input_types, r=r, verbose=False)
name_out_path = os.path.join(outdir, f'gene_names_t{test_fold}.lst')
n = open(name_out_path, 'w')
for ID in gene_names:
    n.write(ID + '\n')
n.close()

save_attribs(input_names=input_names, inputs=inputs, attr=attr, outdir=outdir, method=f't{test_fold}_GradientExplainer', for_modisco=False)


