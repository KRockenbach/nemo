#!/usr/bin/env bash

##################################################################################
#
# MIT License
#
# Copyright (c) 2026 Kevin Rockenbach, Agnieszka Golicz
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
#title: modisco.sh
#description: finds seqlets in 100 bp bins around the TSS and associates them with TF families
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2026-02-22
#version: 1.0.0
#usage: bash modisco.sh
#notes: output is tsv with denoting the TF family and the start and stop coordinate of the respective seqlet,
#       however coordinates are highly biased towards the center of the respective bin
#=========================================================================================================


source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

for ORGANISM in "Bnapus" "Athaliana"; do
  ATTRIB_PATH="../results/nemo/${ORGANISM}/masked_graphpart/attribs"
  mamba activate nemo
  python -m nemo.modisco.transform_arrays $ATTRIB_PATH
  mamba deactivate
  MODISCO_PATH=$(echo $ATTRIB_PATH | sed "s/attribs/modisco/")
  mkdir -p $MODISCO_PATH
  mamba activate modisco
  for SEQ in "promoter" "terminator"; do
    for START in {-500..400..100}; do
      END=$(($START + 100))
      SEQ_FILE="${ATTRIB_PATH}/${SEQ}_seqs_GradientExplainer.modisco_%${START}to${END}%.npz"
      SHAP_FILE="${ATTRIB_PATH}/${SEQ}_shap_GradientExplainer.modisco_%${START}to${END}%.npz"
      ## input sequences are 100 bp intervals
      modisco motifs -w 100 -s $SEQ_FILE -a $SHAP_FILE -n 20000 -o "${MODISCO_PATH}/${SEQ}_%${START}to${END}%_modisco_results.h5"
      modisco report --write-tomtom -i "${MODISCO_PATH}/${SEQ}_%${START}to${END}%_modisco_results.h5" -o "${MODISCO_PATH}/${SEQ}_%${START}to${END}%/" -s "${MODISCO_PATH}/${SEQ}_%${START}to${END}%/" -m ../data/motifs/JASPAR2024_CORE_plants_non-redundant_pfms_meme.txt
    done
  done
  mamba deactivate
done

# associate JASPAR motifs with families
printf "MOTIF\tFAMILY\n" >../data/motifs/motif_families.tsv
for FAMILY in C2C2gata G2like NAC Trihelix bHLH MADS MYB bZIP TCP WRKY C2C2dof Homeobox MYBrelated AP2EREBP HSF; do
  for MOTIF in $(grep "^>" ../data/motifs/${FAMILY}_pfm.txt | cut -f1 -d " " | sed 's/^>//g'); do printf "${MOTIF}\t${FAMILY}\n" >>../data/motifs/motif_families.tsv; done
done
