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
title: shuffle_graphpart_clusters.py
description: shuffles graphpart output for random partitioning
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      python nemo/preprocessing/shuffle_graphpart_clusters.py <graphpart_output>
=========================================================================================================
'''

import sys, os
import pandas as pd

graphpart = pd.read_table(sys.argv[1], index_col=False, header=0, sep=',')
graphpart['cluster'] = graphpart['cluster'].sample(frac=1).values
outpath=sys.argv[1].replace(".csv", "_shuff.csv")
graphpart.to_csv(outpath, sep=',', na_rep='', header=True, index=False)
