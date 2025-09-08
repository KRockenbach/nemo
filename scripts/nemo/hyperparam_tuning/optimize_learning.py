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
title: optimize_learning.py
description: runs fourth step in hyperparameter optimization
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.hyperparam_tuning.optimize_learning <name of study> <number trials before periodic restart> <total number of trials to complete>
notes: run within nemo environment
       tuning restarts periodically to refresh memory
=========================================================================================================
'''

import sys, h5py, os
from math import floor, isnan, ceil
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.optimizers.legacy import Nadam
from tensorflow.keras.losses import Huber
from tensorflow.keras.metrics import R2Score
from tensorflow.keras.callbacks import CSVLogger, EarlyStopping, TerminateOnNaN, LambdaCallback
from tensorflow.keras import backend as K
import optuna
from optuna.integration import TFKerasPruningCallback
from optuna.samplers import TPESampler
from optuna.trial import TrialState
from ..utils.model_utils import *
from ..utils.one_cycle_scheduler_tf.one_cycle_tf.one_cycle_scheduler import OneCycle
from tensorflow_addons.optimizers.weight_decay_optimizers import *


# directory containing training data
datadir = os.path.join("..", "data", "Bnapus", "masked_graphpart_Bn_fold_data")
study_name = sys.argv[1]

# definition and creation of output directory
# outdir with config needs do be created beforehand
outdir = os.path.join("..", "model_configs", study_name)

storage = optuna.storages.RDBStorage(
    url="sqlite:///{}/{}.db".format(outdir, study_name),
    heartbeat_interval=30
)

# define study object
sampler = TPESampler(consider_prior=True,
                     multivariate=True,
                     group=True,
                     n_startup_trials=100,
                     constant_liar=True) # for distributed optimization

study = optuna.create_study(sampler = sampler,
            study_name = study_name,
            direction = "maximize",
            pruner = optuna.pruners.NopPruner(),
            load_if_exists = True,
            storage = storage
)

# CL INPUT 3
n_trials_per_run = int(sys.argv[2])
n_trials_total = int(sys.argv[3]) # should be enough to compensate for startup trials

# define model config get data for the training and validation set depending on config
#########################
input_names=["promoter", "terminator"]

conf_path = os.path.join(outdir, "config.tsv")
config = dict_from_tsv(conf_path)

architecture_path = os.path.join('..', 'model_configs', 'architecture_finetuned', 'hyperparams.tsv')
architecture = dict_from_tsv(architecture_path)

test_fold=0
valid_fold=1


# get full-length (default) data for the training and validation set
# input sequences will be subsampled within the objective function
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

    epochs = trial.suggest_int('epochs', 3, 30, log=True)
    batch_size = 130

    steps_per_epoch = ceil(n_training / batch_size)
    print("Number of steps per epoch = " + str(steps_per_epoch))

    # Nadam_beta_1 set via scheduler
    Nadam_1_minus_beta_2 = trial.suggest_float('Nadam_1_minus_beta_2', 1e-5, 1e-1, log=True)
    Nadam_beta_2 = 1 - Nadam_1_minus_beta_2
    Nadam_epsilon = trial.suggest_float('Nadam_epsilon', 1e-8, 1e-6, log=True)
    Nadam_weight_decay = trial.suggest_float('Nadam_weight_decay', 1e-6, 1e-4, log=True)

    shift_peak = trial.suggest_float('LR_shift_peak', 0.2, 0.8)
    final_lr_scale = trial.suggest_float('LR_final_scale', 0.0, 1.0)
    cycle_size = epochs * steps_per_epoch

    max_lr_lower = 0.01584893 # lowest train_loss overall
    max_lr_upper = 0.37850893 # highest LR before rapid divergence of train_loss
    max_lr = trial.suggest_float("max_lr", max_lr_lower, max_lr_upper, log=True)

    LR_schedule = OneCycle(
        initial_learning_rate=(max_lr/25),
        maximal_learning_rate=max_lr,
        cycle_size=cycle_size,
        shift_peak=shift_peak,
        final_lr_scale=final_lr_scale
    )


    # Nadam_beta_1
    max_momentum = trial.suggest_float('max_momentum', 0.9, 0.99, log=False)
    min_momentum = trial.suggest_float('min_momentum', 0.8, 0.9, log=False)

    momentum_schedule = OneCycle(initial_learning_rate=max_momentum,
                                 maximal_learning_rate=min_momentum,
                                 cycle_size=cycle_size,
                                 shift_peak=shift_peak,
                                 final_lr_scale=1.0
                                 )

    # loss
    huber_delta = trial.suggest_float('huber_delta', 0.135, 13.5)
    loss_fxn = Huber(delta=huber_delta)


    # last layer dropout
    previous_drop = float(architecture["dl_drop_1of2"])
    new_drop = trial.suggest_float("dl_drop_1of2", 0.0, previous_drop)

    # add architecture params
    params = trial.params
    for key, value in architecture.items():
        params[key] = value
    params['batch_size'] = batch_size
    params['dl_drop_1of2'] = new_drop
    print(params)


    model = build_model(params, train, input_names=["promoter", "terminator"])

    # To implement OneCycle policy for LR and momentum, while also using weight decay
    # legacy version of Nadam has to be used, as only this version has ._set_hyper method.
    # However, legacy Nadam does not have weight decay, which is therefore extended using a tensorflow addon
    NadamW = extend_with_decoupled_weight_decay(Nadam)
    optimizer = NadamW(learning_rate=(max_lr/25), # starting value of scheduler
                      beta_1=max_momentum, # starting value of scheduler
                      beta_2=Nadam_beta_2,
                      epsilon=Nadam_epsilon,
                      weight_decay=Nadam_weight_decay)

    optimizer._set_hyper("learning_rate", lambda: LR_schedule(optimizer.iterations))
    optimizer._set_hyper("beta_1", lambda: momentum_schedule(optimizer.iterations))

    model.compile(optimizer=optimizer,
                  loss = loss_fxn,
                  metrics = [R2Score()]
    )


    # custiom callbakc for reporting learning rate and momentum
    class report_LR:
        def __init__(self, model):
            self.model = model
        def on_batch_end(self, batch, logs):
            if batch % 10 == 0:
                logs["lr"] = self.model.optimizer.lr
                logs["beta_1"] = self.model.optimizer.beta_1
    rlr = report_LR(model)
    rlr_cb = LambdaCallback(on_batch_end=lambda batch,
                             logs: rlr.on_batch_end(batch, logs))


    # should be called before training, so that it will be printed even when pruning occurs
    n_pruned_trials = len(study.get_trials(deepcopy=False, states=[TrialState.PRUNED]))
    n_complete_trials = len(study.get_trials(deepcopy=False, states=[TrialState.COMPLETE]))
    n_failed_trials = len(study.get_trials(deepcopy=False, states=[TrialState.FAIL]))
    n_not_failed = n_pruned_trials + n_complete_trials

    if n_not_failed >= n_trials_total:
        exit_function(study)

    print("STUDY STATISTICS" + "\n" + "\n" +
          "Total number of trials:" + "\t" + str(len(study.trials)) + "\n" +
          "Number of pruned trials:" + "\t" + str(n_pruned_trials) + "\n" +
          "Number of failed trials:" + "\t" + str(n_failed_trials) + "\n" +
          "Number of complete trials:" + "\t" + str(n_complete_trials) + "\n")
    if n_complete_trials > 0:
        print("Best score so far:" + "\t" + str(study.best_value) + "\n")


    #define callbacks
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
                 batch_size = batch_size,
                 epochs = epochs,
                 validation_data = ([valid[i] for i in input_names], valid["output"]),
                 callbacks = [csvlog, term, rlr_cb]
    )

    # value used for sampling search space
    score = max(result.history['val_r2_score'])

    return score

study.optimize(objective,
               n_trials = n_trials_per_run, # maximum number of trials
               gc_after_trial = True, # garbage collection after each trial
               catch=(ValueError,  # raised when hyperparameters are incompatible
                      MemoryError) # raised when trial uses too much memory
)








