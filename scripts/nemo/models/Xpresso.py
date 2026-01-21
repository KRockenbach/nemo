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
title: Xpresso.py
description: convenient build function for the Xpresso model
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-01-20
version: 1.0.0
usage:
      (1) run a given script (e.g. script.py) as: python -m script
      (2) within script, import functions using: from ..models.Xpresso import build_nemo
=========================================================================================================
'''

import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, MaxPooling1D, Flatten, Concatenate, Dense, Dropout, LeakyReLU



def build_xpresso():

    ini = "glorot_normal"

    promoter = Input(shape=(10500, 4), name='promoter')

    P = Conv1D(128, 6, padding="same", kernel_initializer=ini)(promoter)
    P = LeakyReLU(alpha=0.1)(P)
    P = MaxPooling1D(30, padding="same")(P)

    P = Conv1D(32, 9, padding="same", kernel_initializer=ini)(P)
    P = LeakyReLU(alpha=0.1)(P)
    P = MaxPooling1D(10, padding="same")(P)

    P = Flatten()(P)
    halflife = Input(shape=(8), name='halflife')
    D = Concatenate()([P, halflife])

    D = Dense(64, kernel_initializer=ini)(D)
    D = LeakyReLU(alpha=0.1)(D)
    D = Dropout(0.00099)(D)

    D = Dense(2, kernel_initializer=ini)(D)
    D = LeakyReLU(alpha=0.1)(D)
    D = Dropout(0.01546)(D)

    D = Dense(1)(D) # output

    return Model(inputs = [promoter, halflife], outputs = D)
