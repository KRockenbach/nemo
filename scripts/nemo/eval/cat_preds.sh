#!/usr/bin/env bash

##################################################################################
#
# MIT License
#
# Copyright (c) 2025 Kevin Rockenbach
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


#=========================================================================================================
#title: cat_preds.sh
#description: concatenate predictions from all ten data folds into one file
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: cat_preds.sh <directory containing predictions>
#=========================================================================================================

PRED_DIR=$1

cat ${PRED_DIR}/actual.t_0.txt > ${PRED_DIR}/concat_actual.lst
cat ${PRED_DIR}/predictions.t_0.txt > ${PRED_DIR}/concat_pred.lst

for TESTFOLD in {1..9}
do
  # only actual*.lst files contain a header
  tail -n+2 ${PRED_DIR}/actual.t_${TESTFOLD}.txt >> ${PRED_DIR}/concat_actual.lst
  tail -n+2  ${PRED_DIR}/predictions.t_${TESTFOLD}.txt >> ${PRED_DIR}/concat_pred.lst
done


echo -e "ID\tActual\tPredicted" > ${PRED_DIR}/concat_preds.tsv
cut -f 2 ${PRED_DIR}/concat_pred.lst | tail -n+2 > ${PRED_DIR}/tmp
tail -n+2 ${PRED_DIR}/concat_actual.lst | paste - ${PRED_DIR}/tmp >> ${PRED_DIR}/concat_preds.tsv
rm  ${PRED_DIR}/tmp
