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
#title: attribution.sh
#description: Calculates attributions for all 10 data folds and concatenates them together for B. napus and A. thaliana
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: bash attribution.sh
#notes:	stores attributions along side the corresponding one-hot sequences for easier importance calculation
#       runs on cpu to parallelize calculation on all 10 folds at once
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo_cpu


for ORGANISM in "Bnapus" "Athaliana"
do
  for TEST in {0..9}
  do
      echo -n "$ORGANISM $TEST "
  done
done | xargs -P 12 -n 2 python3 -m nemo.attribution.attribution

wait

for ORGANISM in "Bnapus" "Athaliana"
do
  ATTR_DIR=../results/nemo/${ORGANISM}/masked_graphpart/attribs
  python3 -m nemo.attribution.concat_attribs $ORGANISM
  cat ${ATTR_DIR}/gene_names_t{0..9}.lst > ${ATTR_DIR}/gene_names.lst
  if [ -f ${ATTR_DIR}/gene_names.lst ]; then
        rm ${ATTR_DIR}/gene_names_t{0..9}.lst
  fi
  for SEQ in "promoter" "terminator"
  do
    for DATA in "shap" "seqs"
    do
      if [ -f ${ATTR_DIR}/${SEQ}_${DATA}_GradientExplainer.full.npz ]; then
        rm ${ATTR_DIR}/${SEQ}_${DATA}_t{0..9}_GradientExplainer.npz
      fi
    done
  done
done

mamba deactivate
