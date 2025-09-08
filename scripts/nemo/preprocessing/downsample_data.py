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
title: downsample_data.py
description: downsamples B. napus data folds to size of corresponding A. thaliana data folds
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      from within nemo/scripts/ directory:
      python nemo/preprocessing/downsample_data.py <Bnapus_data_directoy> <Output_data_directoy>
=========================================================================================================
'''

import os, sys
import pandas as pd

datadir=sys.argv[1]
outdir=sys.argv[2]
Ath_GP_lst = datadir.split("/")[0:-2]
Ath_GP_lst.extend(["Athaliana", "derived_data", "graphpart_result.csv"])
Ath_GP_path = os.path.join(*Ath_GP_lst)
GP_table = pd.read_csv(Ath_GP_path, delimiter=",", header=0)

for fold in range(10):

    path = os.path.join(datadir, f'fold_{fold}.feather')
    fold_table = pd.read_feather(path)
    N = sum(GP_table.loc[:,"cluster"]==fold)
    ds_table = fold_table.sample(n=N, axis=0)
    ds_table.to_feather(os.path.join(outdir, f"fold_{fold}.feather"))
