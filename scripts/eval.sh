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
#title: eval.sh
#description: Creates prediction sets from trained models, calculates and aggregates performances and classifies genes based on expression and prediction
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: bash eval.sh
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo

ROOT=".."

echo -e "model\ttest_organism\tmasking\tpartitioning\ttrain_organism\ttest_fold\tvalid_fold\tN\trsq" > "${ROOT}/results/rsq_df.tsv"


# perfomance evaluation
# models parallelized
ORGANISM="Bnapus"
MASKING="masked"
PARTITIONING="graphpart_Bn"
FOLDDIR=${ROOT}"/data/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}"_fold_data"
NEMO_WEIGHTDIR=${ROOT}"/model_weights/nemo/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}
XPRESSO_WEIGHTDIR=${ROOT}"/model_weights/xpresso/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}
CUDA_VISIBLE_DEVICES=0 python3 -m nemo.eval.eval "$FOLDDIR" "$NEMO_WEIGHTDIR" & pid1=$!
CUDA_VISIBLE_DEVICES=1 python3 -m nemo.eval.eval "$FOLDDIR" "$XPRESSO_WEIGHTDIR" & pid2=$!
wait $pid1 $pid2

for MODEL in "nemo" "xpresso"
do
  ORGANISM="Bnapus"
  MASKING="masked"
  PARTITIONING="graphpart_Bn"
  #                                  <model_name> <test_organism> <test_masking> <partitioning> <train_organism> <train_masking>
  Rscript nemo/eval/get_rsq_matrix.R $MODEL $ORGANISM $MASKING $PARTITIONING $ORGANISM
done


## performance comparisons w.r.t. masking and paritioning (only Bnapus)
# partitioning parallelized
for MASKING in "clear" "masked" "only_cds"
do
  GP_TESTDIR=${ROOT}"/data/Bnapus/"${MASKING}"_graphpart_fold_data"
  RD_TESTDIR=${ROOT}"/data/Bnapus/"${MASKING}"_random_fold_data"
  GP_WEIGHTDIR=${ROOT}"/model_weights/nemo/Bnapus/"${MASKING}"_graphpart"
  RD_WEIGHTDIR=${ROOT}"/model_weights/nemo/Bnapus/"${MASKING}"_random"
  CUDA_VISIBLE_DEVICES=0 python3 -m nemo.eval.eval "$GP_TESTDIR" "$GP_WEIGHTDIR" & pid1=$!
  CUDA_VISIBLE_DEVICES=1 python3 -m nemo.eval.eval "$RD_TESTDIR" "$RD_WEIGHTDIR" & pid2=$!
  wait $pid1 $pid2
done


## cross-species comparisons
# test organism parallelized
for MODELORG in "Bnapus" "Athaliana"
do
  # Bnapus test set downsampled to size of Athaliana test sets
  Bn_TESTDIR=${ROOT}"/data/Bnapus/masked_graphpartDS_fold_data"
  At_TESTDIR=${ROOT}"/data/Athaliana/masked_graphpart_fold_data"
  WEIGHTDIR=${ROOT}"/model_weights/nemo/"${MODELORG}"/masked_graphpart"
  CUDA_VISIBLE_DEVICES=0 python3 -m nemo.eval.eval "$Bn_TESTDIR" "$WEIGHTDIR" & pid1=$!
  CUDA_VISIBLE_DEVICES=1 python3 -m nemo.eval.eval "$At_TESTDIR" "$WEIGHTDIR" & pid2=$!
  wait $pid1 $pid2

  for TESTORG in "Bnapus" "Athaliana"
  do
    if [[ $TESTORG == "Bnapus" && $MODELORG == "Bnapus" ]]; then
      for MASKING in "clear" "masked" "only_cds"
      do
        for PARTITIONING in "graphpart" "random"
        do
          #                                  <model_name> <test_organism> <masking> <partitioning> <train_organism>
          Rscript nemo/eval/get_rsq_matrix.R "nemo" "Bnapus" $MASKING $PARTITIONING "Bnapus"
        done
      done
    else
      #                                  <model_name> <test_organism> <masking> <partitioning> <train_organism>
      Rscript nemo/eval/get_rsq_matrix.R "nemo" $TESTORG "masked" "graphpart" $MODELORG
    fi
  done
done

## classify based on expression, specificity and prediction
for ORGANISM in "Bnapus" "Athaliana"
do
  EXPRESSION=${ROOT}"/data/"${ORGANISM}"/derived_data/expr_matrix.tsv"
  MODEL="nemo"
  PARTITIONING="graphpart"
  MASKING="masked"
  INITIAL=$(echo $ORGANISM | head -c 1)
  PREDDIR=${ROOT}"/results/"${MODEL}"/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}"/"${MODEL}${INITIAL}"_preds"

  bash nemo/eval/cat_preds.sh $PREDDIR
  EXPR_IDDIR=${ROOT}"/results/"${MODEL}"/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}"/IDs/expression"
  SPEC_IDDIR=${ROOT}"/results/"${MODEL}"/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}"/IDs/specificity"
  PRED_IDDIR=${ROOT}"/results/"${MODEL}"/"${ORGANISM}"/"${MASKING}"_"${PARTITIONING}"/IDs/prediction"
  Rscript nemo/eval/over_under_pred.R $PREDDIR $PRED_IDDIR
  Rscript nemo/eval/get_expression_class_lists.R $PREDDIR $EXPR_IDDIR
  Rscript nemo/eval/get_specificity_class_lists.R $EXPRESSION $SPEC_IDDIR
done

mamba deactivate
