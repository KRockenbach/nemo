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
#date: 2026-01-21
#version: 1.0.1
#usage: bash train_Bnapus.sh
#notes: To train in parallel device 1 is used for A. thaliana, device 0 used for B. napus
#=========================================================================================================

DEVICE=0

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo

for MASK in "masked" "clear" "only_cds"
do
  for PARTITION in "graphpart" "random"
  do
    ORGANISM="Bnapus"
    TPM_TYPE="median"
    for TEST in {0..9}
    do
      VALID="None"
      N="None"
      CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
    done
    if [[ $MASK == "masked" && $PARTITION == "graphpart" ]]; then
      for N in {0..9}
      do
        for TPM_TYPE in "max" "median"
        do
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_Xpresso $N $TPM_TYPE
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_Basenji-5K $N $TPM_TYPE
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_PhytoExpr_CNN $N $TPM_TYPE
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_PhytoExpr_transformer $N $TPM_TYPE
          TEST=0
          VALID=1
          PARTITION="graphpart_Bn"
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
          VALID="None"
          CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
        done
      done
      # train on full set
      VALID="None"
      TEST="None"
      N="None"
      PARTITION="graphpart"
      TPM_TYPE="median"
      CUDA_VISIBLE_DEVICES=$DEVICE python -m nemo.training.train_nemo $ORGANISM $MASK $PARTITION $TEST $VALID $N $TPM_TYPE
    fi
  done
done

mamba deactivate
