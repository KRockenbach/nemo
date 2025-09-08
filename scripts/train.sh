#!/usr/bin/env bash

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


#========================================================================================================
#title: run_setup.sh
#description: trains model for a given species, masking-method and partitioning method, while holding out defined test/validation sets
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage:
#      CUDA_DEVICE=0
#      for MASK in masked clear
#      do
#        for PARTITION in graphpart random
#        do
#          for TEST in {0..9}
#          do
#            for VALID in {0..9}
#            do
#              for N in {0..9}
#              do
#                bash train.sh $CUDA_DEVICE Bnapus nemo $MASK $PARTITION $TEST $VALID $N
#              done
#            done
#          done
#        done
#      done
#
#notes: to run training on full set, set TEST_SET and VALID_SET both to "None"
#       to run training on 90% set, set only VALID_SET to "None"
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

export CUDA_VISIBLE_DEVICES=$1

ORGANISM=$2
MODEL=$3
MASK=$4
PARTITION=$5
TEST_SET=$6
VALID_SET=$7
N=$8


mamba activate nemo
python nemo/run_training.py $ORGANISM $MODEL $MASK $PARTITION $TEST_SET $VALID_SET $N
mamba deactivate

