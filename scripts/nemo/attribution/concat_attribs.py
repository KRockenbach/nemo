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
title: concat_attribs.py
description: concatenates the attributions/sequences of all test folds into one file
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.attribution.concat_attribs <species>
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


root = '..'

organism = sys.argv[1] # directory containing training data
model_name = "nemo"
modeldir = os.path.join(root, 'model_configs', model_name)
resultdir = modeldir.replace('model_configs', 'results')
conf_path = os.path.join(modeldir, 'config.tsv')
config = dict_from_tsv(conf_path)
outdir = os.path.join(resultdir, organism, 'masked_graphpart', 'attribs')


##########################################################

input_names = []
if to_bool(config['use_promoter']):
    input_names.append('promoter')
if to_bool(config['use_terminator']):
    input_names.append('terminator')
if to_bool(config['use_halflife']):
    input_names.append('halflife')


attr = []
inputs = []
for in_name in input_names:
    for test_fold in range(10):
        attr_path = os.path.join(outdir, in_name + f"_shap_t{test_fold}_GradientExplainer.npz")
        in_path = os.path.join(outdir, in_name + f"_seqs_t{test_fold}_GradientExplainer.npz")
        new_attr_arr = np.load(attr_path)['arr_0']
        new_input_arr = np.load(in_path)['arr_0']
        if test_fold == 0:
            attr_arr = new_attr_arr
            input_arr = new_input_arr
        else:
            attr_arr = np.concatenate((attr_arr, new_attr_arr), axis=0)
            input_arr = np.concatenate((input_arr, new_input_arr), axis=0)
    attr.append(attr_arr)
    inputs.append(input_arr)
save_attribs(input_names=input_names, inputs=inputs, attr=attr, outdir=outdir, method=f'GradientExplainer.full', for_modisco=False)



