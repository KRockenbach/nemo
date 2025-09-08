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
#date: 2025-08-28
#version: 1.0.0
#usage: bash train_Bnapus.sh
#notes: To train in parallel device 1 is used for A. thaliana, device 0 used for B. napus
#=========================================================================================================

DEVICE=0

for MASK in "masked" "clear" "only_cds"
do
  for PARTITION in "graphpart" "random"
  do
    ORGANISM="Bnapus"
    MODEL="nemo"
    for TEST in {0..9}
    do
      VALID="None"
      N="None"
      bash train.sh $DEVICE $ORGANISM $MODEL $MASK $PARTITION $TEST $VALID $N
    done
    if [[ $MASK == "masked" && $PARTITION == "graphpart" ]]; then
      MODEL="xpresso"
      for N in {0..9}
      do
        TEST=0
        VALID=1
        bash train.sh $DEVICE $ORGANISM $MODEL $MASK "graphpart_Bn" $TEST $VALID $N
      done
      MODEL="nemo"
      for N in {0..9}
      do
        TEST=0
        VALID=1
        bash train.sh $DEVICE $ORGANISM $MODEL $MASK "graphpart_Bn" $TEST $VALID $N
        TEST=0
        VALID="None"
        bash train.sh $DEVICE $ORGANISM $MODEL $MASK "graphpart_Bn" $TEST $VALID $N
      done
      # train on full set
      VALID="None"
      TEST="None"
      N="None"
      bash train.sh $DEVICE $ORGANISM "nemo" "masked" "graphpart" $TEST $VALID $N
    fi
  done
done
