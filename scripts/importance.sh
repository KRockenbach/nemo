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
#title: importance.sh
#description: Calculates group-specific importance (trimmed = +- 1kb around TSS/TTS; full = entire input)
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: bash importance.sh
#notes: group directories contain lists of gene IDs for the respective groups
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo

MODEL="nemo"
for ORGANISM in "Athaliana" "Bnapus"
do
  GROUP_DIR="../data/${ORGANISM}/parent_data/TF_ids/"
  CLUSTERING="TF"
  python3 -m nemo.attribution.get_trimmed_importance $GROUP_DIR $CLUSTERING $MODEL $ORGANISM

  GROUP_DIR="../results/${MODEL}/${ORGANISM}/masked_graphpart/IDs/expression/"
  CLUSTERING="expression"
  python3 -m nemo.attribution.get_full_importance $GROUP_DIR $CLUSTERING $MODEL $ORGANISM

  GROUP_DIR="../results/${MODEL}/${ORGANISM}/masked_graphpart/IDs/prediction/"
  CLUSTERING="prediction"
  python3 -m nemo.attribution.get_full_importance $GROUP_DIR $CLUSTERING $MODEL $ORGANISM

done

mamba deactivate
