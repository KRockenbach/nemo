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
title: nemo.py
description: convenient build function for the nemo model
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-01-20
version: 1.0.0
usage:
      (1) run a given script (e.g. script.py) as: python -m script
      (2) within script, import functions using: from ..models.nemo import build_nemo
=========================================================================================================
'''

import os, sys
from math import ceil, floor
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, AveragePooling1D, Flatten, Concatenate, Dense, BatchNormalization, Dropout
from tensorflow_addons.layers import GELU


def build_nemo(promoter, temrinator):

    ini = "glorot_normal"
    def get_stride(fraction, size):
        return max([1, ceil(fraction*size)])

    # promoter branch
    P = Conv1D(288, 5, padding = 'same', kernel_initializer = ini)(promoter)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(106, 8, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)
    P = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(P)

    P = Conv1D(227, 9, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(187, 10, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)

    P = Conv1D(64, 46, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)
    P = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(P)

    P = Conv1D(152, 62, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=0.81669)(P)
    P = GELU(approximate=False)(P)
    P = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(P)

    P = Flatten()(P)

    # terminator branch
    T = Conv1D(219, 4, padding = 'same', kernel_initializer = ini)(terminator)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(449, 4, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)
    T = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(T)

    T = Conv1D(259, 16, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(154, 15, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)

    T = Conv1D(105, 20, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)
    T = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(T)

    T = Conv1D(152, 22, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=0.81669)(T)
    T = GELU(approximate=False)(T)
    T = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(T)

    T = Flatten()(T)

    # dense layers
    D = Concatenante(axis=1)([P,T])
    D = Dense(750, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)

    D = Dense(3, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=0.81669)(D)
    D = GELU(approximate=False)(D)
    D = Dropout(0.00124)(D)
    D = Dense(1)(D) # single-regression output
    return Model(inputs = [promoter, terminator], outputs = D)
