#!/usr/bin/python3

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


'''
=========================================================================================================
title: validate_architecture.py
description: runs second step in hyperparameter optimization
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.hyperparam_tuning.validate_architecture <name of study>
notes: run within nemo environment
       tuning restarts periodically to refresh memory
=========================================================================================================
'''

import sys, h5py, os
from math import floor, isnan, ceil
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers.experimental import Adam
from tensorflow.keras.losses import Huber
from tensorflow.keras.metrics import R2Score
from tensorflow.keras.optimizers.schedules import CosineDecayRestarts
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, TerminateOnNaN, Callback
from tensorflow.keras import backend as K
import optuna
from optuna.integration import TFKerasPruningCallback
from optuna.samplers import TPESampler
from optuna.trial import TrialState
from ..utils.model_utils import *


# CL INPUT 1
# directory containing training data
organism = "Bnapus"
datadir = os.path.join("..", "data", organism, "masked_graphpart_Bn_fold_data")

# CL INPUT 2
study_name = sys.argv[1]

# definition and creation of output directory
# CL INPUT 1
# outdir with config needs do be created beforehand
outdir = os.path.join("..", "model_configs", study_name)

os.makedirs(outdir, exist_ok=True)

storage = optuna.storages.RDBStorage(
    url="sqlite:///{}/{}.db".format(outdir, study_name),
    heartbeat_interval=30
)

# define study object
n_epo = 10 # max number of epochs
pruner = optuna.pruners.NopPruner()
search_space = {
    'P_branch': ['P_params', 'T_params'],
    'T_branch' : ['P_params', 'T_params'],
    'default_conv_stride' : [True, False],
    'default_pool_stride' : [True, False],
    'pooling': ['optimized', 'avrg', 'max']
} # 96 combinations

sampler = optuna.samplers.GridSampler(search_space)

study = optuna.create_study(sampler = sampler,
            study_name = study_name,
            direction = "minimize", # val_loss will be maximized
            pruner = pruner,
            load_if_exists = True,
            storage = storage
)

# define model config get data for the training and validation set depending on config
#########################
input_names=["promoter", "terminator"]

conf_path = os.path.join(outdir, "config.tsv")
config = dict_from_tsv(conf_path)

test_fold=0
valid_fold=1


# input sequences preloaded and subsampled to required length
#########################
train = get_set(config, outP=5000, outT=5000, inP=1200, inT=1200, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
valid = get_set(config, outP=5000, outT=5000, inP=1200, inT=1200, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)

# get number of training examples
n_training = train[input_names[0]].shape[0]
print("Number of training examples = " + str(n_training))

###############################################################################
###############################################################################




def exit_function(study):
    # log results of optimization
    # save to hyperparameter output file
    ###########################################################################################
    df = study.trials_dataframe(attrs=('number', 'value', 'state', 'params'), multi_index=False)
    df.to_csv(os.path.join(outdir, 'trial_results.csv'), index=False, header=True, sep='\t')


    n_pruned_trials = len(study.get_trials(deepcopy=False, states=[TrialState.PRUNED]))
    n_complete_trials = len(study.get_trials(deepcopy=False, states=[TrialState.COMPLETE]))
    n_failed_trials = len(study.get_trials(deepcopy=False, states=[TrialState.FAIL]))
    bestTrial = study.best_trial

    summary = open(os.path.join(outdir, 'optimization_summary.txt'), 'w')
    summary.write("STUDY STATISTICS" + "\n" + "\n" +
                  "Total number of trials:" + "\t" + str(len(study.trials)) + "\n" +
                  "Number of pruned trials:" + "\t" + str(n_pruned_trials) + "\n" +
                  "Number of complete trials:" + "\t" + str(n_complete_trials) + "\n" +
                  "Number of failed trials:" + "\t" + str(n_failed_trials) + "\n" +
                  "Best trial:" + "\n" + "\n" + str(study.best_trial) + "\n" +
                  "Value:" + "\t" + str(bestTrial.value)+"\n" + "\n" +
                  "HYPERPARAMETERS:" + "\n" + "\n")

    best_params = bestTrial.params

    for key, value in best_params.items():
        summary.write("\t" + "{}:".format(key) + "\t" + "{}".format(value) + "\n")
    summary.close()

    bestHPs = open(os.path.join(outdir, 'hyperparams.tsv'), 'w')
    for key, value in bestTrial.params.items():
        bestHPs.write("{}".format(key)+"\t" + "{}".format(value) + "\n")
    bestHPs.close()

    print("Run Successful")
    print("STUDY STATISTICS" + "\n" + "\n" +
          "Total number of trials:" + "\t" + str(len(study.trials)) + "\n" +
          "Number of pruned trials:" + "\t" + str(n_pruned_trials) + "\n" +
          "Number of failed trials:" + "\t" + str(n_failed_trials) + "\n" +
          "Number of complete trials:" + "\t" + str(n_complete_trials) + "\n")
    exit(0)




##############################################################################
##############################################################################


#Function to build model with trial-specific hyperparameters
def objective(trial):

    K.clear_session() # refresh memory for each trial

    ###############################
    # Define hyperparameters
    ###############################

    trial.suggest_categorical('P_branch', ['P_params', 'T_params'])
    trial.suggest_categorical('T_branch', ['P_params', 'T_params'])
    trial.suggest_categorical('default_conv_stride', [True, False])
    trial.suggest_categorical('default_pool_stride', [True, False])
    trial.suggest_categorical('pooling', ['optimized', 'avrg', 'max'])

    fixed_params = {
        'outside_P': 5000,
        'outside_T': 5000,
        'inside_P': 1200,
        'inside_T': 1200,
        'batch_size': 32,
        'use_drop': True,
        'use_conv_drop': False,
        'use_lstm_drop': False,
        'use_exp': False,
        'ndl': 1,
        'dl_dim_0of1': 8,
        'dl_drop_0of1': 0.06766238081567807,
        'activation': 'gelu',
        'kernel_ini_gelu': 'glorot_normal',
        'use_Bnorm': True,
        'Bnorm_momentum': 0.9176885178453292,
        'skip_type': 'none'
    }

    branch_P_params_P = {
        'ncl_P': 6,
        'nlstm_P': 0,
        'filter_P_0of6': 512, # 4**4 = 256 --> twice as many filters as 4-mers
        'kernel_P_0of6': 4,
        'pool_size_P_0of6': 4,
        'conv_frac_P_0of6': 0.0,
        'pool_type_P_0of6': 'max',
        'pool_frac_P_0of6': 0.7,
        'filter_P_1of6': 256,
        'kernel_P_1of6': 4,
        'pool_size_P_1of6': 64,
        'conv_frac_P_1of6': 0.2,
        'pool_type_P_1of6': 'avrg',
        'pool_frac_P_1of6': 0.5,
        'filter_P_2of6': 128,
        'kernel_P_2of6': 16,
        'pool_size_P_2of6': 2,
        'conv_frac_P_2of6': 0.2,
        'pool_type_P_2of6': 'avrg',
        'pool_frac_P_2of6': 1.0,
        'filter_P_3of6': 128,
        'kernel_P_3of6': 4,
        'pool_size_P_3of6': 2,
        'conv_frac_P_3of6': 0.5,
        'pool_type_P_3of6': 'max',
        'pool_frac_P_3of6': 0.9,
        'filter_P_4of6': 512,
        'kernel_P_4of6': 32,
        'pool_size_P_4of6': 8,
        'conv_frac_P_4of6': 0.5,
        'pool_type_P_4of6': 'avrg',
        'pool_frac_P_4of6': 0.5,
        'filter_P_5of6': 64,
        'kernel_P_5of6': 32,
        'pool_size_P_5of6': 4,
        'conv_frac_P_5of6': 0.3,
        'pool_type_P_5of6': 'max',
        'pool_frac_P_5of6': 0.7
    }

    branch_T_params_P = {}
    for key in branch_P_params_P:
        new_key = key.replace('_P','_T')
        branch_T_params_P[new_key]=branch_P_params_P[key]

    branch_T_params_T = {
        'ncl_T': 2,
        'nlstm_T': 2,
        'filter_T_0of2': 32,
        'kernel_T_0of2': 4,
        'pool_size_T_0of2': 64,
        'conv_frac_T_0of2': 0.5,
        'pool_type_T_0of2': 'avrg',
        'pool_frac_T_0of2': 0.7,
        'filter_T_1of2': 32,
        'kernel_T_1of2': 8,
        'pool_size_T_1of2': 4,
        'conv_frac_T_1of2': 0.0,
        'pool_type_T_1of2': 'max',
        'pool_frac_T_1of2': 1.0,
        'lstm_units_T_0of2': 32,
        'lstm_units_T_1of2': 128
    }

    branch_P_params_T = {}
    for key in branch_T_params_T:
        new_key = key.replace('_T','_P')
        branch_P_params_T[new_key]=branch_T_params_T[key]

    # merge dictionaries
    params = fixed_params
    print(trial.params)
    if trial.params['P_branch'] == 'P_params':
        params = params | branch_P_params_P
    else:
        params = params | branch_P_params_T
    if trial.params['T_branch'] == 'T_params':
        params = params | branch_T_params_T
    else:
        params = params | branch_T_params_P
    if trial.params['default_conv_stride']:
        for key in params:
            if key.startswith('conv_frac'):
                params[key] = 0.0
    if trial.params['default_pool_stride']:
        for key in params:
            if key.startswith('pool_frac'):
                params[key] = 1.0
    if trial.params['pooling'] == 'max':
        for key in params:
            if key.startswith('pool_type'):
                params[key] = 'max'
    elif trial.params['pooling'] == 'avrg':
        for key in params:
            if key.startswith('pool_type'):
                params[key] = 'avrg'

    ###################################################################################
    # end of search space definition
    ###################################################################################

    ###############################################################################
    # define model
    ###############################################################################

    print(params)

    model = build_model(params, train, input_names=["promoter", "terminator"])

    ###############################################################################
    # end of model definition
    ###############################################################################


    # same optimizer settings as in Xpresso
    loss_fxn = "mse"
    optimizer_fxn = Adam(learning_rate=0.001,
                             beta_1=0.9,
                             beta_2=0.999,
                             epsilon=1e-08,
                             weight_decay=0.0)

    model.compile(optimizer_fxn,
                        loss = loss_fxn,
                        metrics = [R2Score()]
    )


    n_pruned_trials = len(study.get_trials(deepcopy=False, states=[TrialState.PRUNED]))
    n_complete_trials = len(study.get_trials(deepcopy=False, states=[TrialState.COMPLETE]))
    n_failed_trials = len(study.get_trials(deepcopy=False, states=[TrialState.FAIL]))


    print("STUDY STATISTICS" + "\n" + "\n" +
          "Total number of trials:" + "\t" + str(len(study.trials)) + "\n" +
          "Number of pruned trials:" + "\t" + str(n_pruned_trials) + "\n" +
          "Number of failed trials:" + "\t" + str(n_failed_trials) + "\n" +
          "Number of complete trials:" + "\t" + str(n_complete_trials) + "\n")
    if n_complete_trials > 0:
        print("Best score so far:" + "\t" + str(study.best_value) + "\n")


    # define callbacks
    # pruning based on intermediate validation loss values (minimization)
    pruning = TFKerasPruningCallback(trial, 'val_loss')


    #Thanks to https://stackoverflow.com/a/37296168
    class StopOnNegativeR2(Callback):
        def __init__(self, monitor='r2_score', value=0.0, verbose=1):
            super(Callback, self).__init__()
            self.monitor = monitor
            self.value = value
            self.verbose = verbose

        def on_epoch_end(self, epoch, logs={}):
            current = logs.get(self.monitor)
            if current is None:
                warnings.warn("Early stopping requires %s available!" % self.monitor, RuntimeWarning)

            if current < self.value:
                if self.verbose > 0:
                    print("Epoch %05d: stopping on negative R2" % epoch)
                self.model.stop_training = True

    negR2 = StopOnNegativeR2()



    try:
        df = study.trials_dataframe()
        n = df.shape[0]  # n_trials
        del df # free memory back up
    except:
        n = 0
        pass
    logdir = os.path.join(outdir, "csv_logs")
    if not os.path.exists(logdir):
        os.makedirs(logdir)
    logfile = os.path.join(logdir,f'full_training_log_{n}.csv')
    csvlog = CSVLogger(logfile, append=True, separator='\t')
    term = TerminateOnNaN()

    # fit model
    result = model.fit(
                 [train[i] for i in input_names],
                 train["output"],
                 batch_size = params['batch_size'],
                 epochs = n_epo,
                 validation_data = ([valid[i] for i in input_names], valid["output"]),
                 callbacks = [pruning, csvlog, term, negR2]
    )

    # value used for sampling search space
    score = min(result.history['val_loss'])

    optlog = open(os.path.join(outdir, 'optimization_log.txt'), 'a')
    for key, value in trial.params.items():
        optlog.write("+++++++++++++++++++++++++")
        optlog.write("\t" + "{}:".format(key) + "\t" + "{}".format(value) + "\n")
        optlog.write("\n")
        optlog.write("Score: " + str(score) + "\n\n")
    optlog.close()

    return score

study.optimize(objective)

exit_function(study)






