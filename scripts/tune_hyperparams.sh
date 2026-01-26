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
#title: tune_hyperparams.sh
#description: wrapper for running hyperparameter optimization studies
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage:
#      bash tune_hyperparams.sh 0 optimize_architecture.py architecture_optimized 1100 200 False
#      CUDA_VISIBLE_DEVICES=0 python -m nemo.hyperparam_tuning.validate_architecture.py architectire_validated
#      bash tune_hyperparams.sh 0 finetune_architecture.py architecture_finetuned 1100 200 False
#      bash tune_hyperparams.sh 0 optimize_learning.py learning_optimized 510 300 True
#      bash tune_hyperparams.sh 0 finetune_architecture.py learning_finetuned 510 300 True
#      bash tune_hyperparams.sh 0 final_touches.py nemo 310 310 True
#      CUDA_VISIBLE_DEVICES=1 python -m nemo.hyperparam_tuning.validate_final_architecture final_architectire_validated
#      bash tune_hyperparams.sh 0 expression_dynamics.py nemo2_foundation 1100 50 True
#      bash tune_hyperparams.sh 0 expression_dynamics_2.py nemo2_finetuning 250 50 True
#notes:
#      Script can be executed multiple times in parallel, if multiple GPUs are available. The first argument should be changed accordingly
#      Parallel optimization will be performed using a common database
#      Depending on search space, failure rate might be high. Script is restarted after N_TRAILS_PER_RUN trials, until TRIALS_TO_COMPLETE is reached
#      Study can be interrupted and resumed any time
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

export CUDA_VISIBLE_DEVICES=$1

SCRIPTNAME=$2
TRUNCATED_SCRIPTNAME=$(echo $SCRIPTNAME | sed 's/.py$//')
SCRIPTNAME=$TRUNCATED_SCRIPTNAME
STUDY_NAME=$3


OUTDIR="../model_configs/"$STUDY_NAME
mkdir -p $OUTDIR


mamba activate nemo

N_TRIALS_COMPLETE=0
TRIALS_TO_COMPLETE=$4 # number of trials for bayesian optimization + number of warmup trials (random sampling)
N_TRIALS_PER_RUN=$5 # interval after which script is restarted, mitigates risk of failure due to memory leakage


# set config
cfg=$OUTDIR"/config.tsv"
echo -e "use_promoter\tTrue" > $cfg
echo -e "use_terminator\tTrue" >> $cfg
echo -e "use_halflife\tFalse" >> $cfg
# True or False
ONE_CYCLE=$6
if [[ $ONE_CYCLE == "True" ]]; then
  echo -e "use_xpresso_optimization\tFalse" >> $cfg
  echo -e "use_one_cycle_optimization\tTrue" >> $cfg
else
  echo -e "use_xpresso_optimization\tTrue" >> $cfg
  echo -e "use_one_cycle_optimization\tFalse" >> $cfg
fi

while [ $N_TRIALS_COMPLETE -lt $TRIALS_TO_COMPLETE ]
do
  python3 -m nemo.hyperparam_tuning.$TRUNCATED_SCRIPTNAME $STUDY_NAME $N_TRIALS_PER_RUN $TRIALS_TO_COMPLETE
  N_TRIALS_COMPLETE=$(python nemo/hyperparam_tuning/get_completed_trials.py $STUDY_NAME)
  echo $N_TRIALS_COMPLETE" trials completed or pruned"
  sleep 5
done

mamba deactivate
