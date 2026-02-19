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
#title: train_Bnapus.sh
#description: runs required training for B. napus
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2026-01-27
#version: 1.0.2
#usage: bash train_Bnapus.sh
#notes: To train in parallel device 1 is used for A. thaliana, device 0 used for B. napus
#=========================================================================================================

DEVICE=1

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo

#for MASK in "masked" "clear" "only_cds"
#do
#  for PARTITION in "graphpart" "random"
#  do
#    N="None"
#    ORGANISM="Bnapus"
#    TPM_TYPE="median"
#    for TEST in {0..9}
#    do
#      VALID="None"
#      CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
#    done
#    # train on full set
#    VALID="None"
#    TEST="None"
#    PARTITION="graphpart"
#    CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
#  done
#done

MASK="masked"
PARTITION="graphpart"
ORGANISM="Bnapus"
PARTITION="graphpart_Bn"
for N in {0..9}
do
  for TPM_TYPE in "max" #"median"
  do
    TEST=0
    VALID=1
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_Xpresso $N $TPM_TYPE
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_Xpresso_no_halflife $N $TPM_TYPE
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_Basenji-5K $N $TPM_TYPE
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_PhytoExpr_CNN $N $TPM_TYPE
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_PhytoExpr_transformer $N $TPM_TYPE
    CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
    #VALID="None"
    #CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
  done
done

mamba deactivate
