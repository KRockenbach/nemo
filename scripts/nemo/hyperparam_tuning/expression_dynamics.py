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
title: final_touches.py
description: runs final step in hyperparameter optimization
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-01-07
version: 1.0.0
usage: python -m nemo.hyperparam_tuning.expression_dynamics <name of study> <number trials before periodic restart> <total number of trials to complete>
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
from ..utils.ModelTools import *
from ..utils.one_cycle_scheduler_tf.one_cycle_tf.one_cycle_scheduler import OneCycle
from tensorflow_addons.optimizers.weight_decay_optimizers import *

# directory containing training data
datadir = os.path.join("..", "data", "Bnapus", "masked_graphpart_Bn_fold_data")

# CL INPUT 1
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
#                     seed=1234,
                     multivariate=True,
                     group=True,
                     n_startup_trials=10,
                     constant_liar=True) # for distributed optimization

study = optuna.create_study(sampler = sampler,
            study_name = study_name,
            direction = "minimize",
            pruner = optuna.pruners.NopPruner(), # no pruning
            load_if_exists = True,
            storage = storage
)

# CL INPUT 3
n_trials_per_run = int(sys.argv[2])
n_trials_total = int(sys.argv[3]) # should be enough to compensate for startup trials + at least 1000

# define model config get data for the training and validation set depending on config
#########################
input_names=["promoter", "terminator"]

conf_path = os.path.join(outdir, "config.tsv")
config = dict_from_tsv(conf_path)

## load nemo hyperparameters
architecture_path = os.path.join('..', 'model_configs', 'nemo', 'hyperparams.tsv')
architecture = dict_from_tsv(architecture_path)
params = {}
for key, value in architecture.items():
    params[key] = value

test_fold=0
valid_fold=1


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
        params[key] = value
        summary.write("\t" + "{}:".format(key) + "\t" + "{}".format(value) + "\n")
    summary.close()

    HP = open(os.path.join(outdir, 'hyperparams.tsv'), 'w')
    for key, value in params:
        HP.write("{}".format(key)+"\t" + "{}".format(value) + "\n")
    HP.close()

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

    # this round of optimization is primarily there to push the target LR as high as possible

    K.clear_session() # refresh memory for each trial

    epochs = 6
    batch_size = 130

    steps_per_epoch = ceil(n_training / batch_size)
    print("Number of steps per epoch = " + str(steps_per_epoch))

    shift_peak = 0.7953957668000962
    final_lr_scale = 0.8757807116665555
    cycle_size = epochs * steps_per_epoch

    max_lr = 0.0663012475147729
    min_lr = 0.0012994663859360526

    Bnorm_momentum = 0.816691400963907
    huber_delta = 2.9465947350775075

    # Nadam_beta_1 set via scheduler
    Nadam_1_minus_beta_2 = 0.0017932649551730192
    Nadam_beta_2 = 1 - Nadam_1_minus_beta_2
    lr_over_epsilon = 0.04872829743137629 / 8.327220255431336e-07
    Nadam_epsilon = max_lr / lr_over_epsilon
    Nadam_weight_decay = 3.209356965523218e-05

    LR_schedule = OneCycle(
        initial_learning_rate=min_lr,
        maximal_learning_rate=max_lr,
        cycle_size=cycle_size,
        shift_peak=shift_peak,
        final_lr_scale=final_lr_scale
    )


    # Nadam_beta_1
    max_momentum = 0.927359606261753
    min_momentum = 0.8899148069769187

    momentum_schedule = OneCycle(initial_learning_rate=max_momentum,
                                 maximal_learning_rate=min_momentum,
                                 cycle_size=cycle_size,
                                 shift_peak=shift_peak,
                                 final_lr_scale=1.0
                                 )

    # loss
    #huber_delta = trial.suggest_float('huber_delta', 1.0, 3.0, log=False) # linear in last round
    loss_fxn = Huber(delta=huber_delta)


    # last layer dropout
    drop = 0.0012428917539097867 # dropout from last round


    dense_dim = trial.suggest_integer("dl_dim_1of2", 5, 30, log=True)
    dropout = trial.suggest_float('dl_drop_1of2', 5e-3, 5e-1, log=True)

    HP_dir = os.path.join(outdir, 'hyperparams')
    os.makedirs(HP_dir, exist_ok=True)
    HP_out = open(os.path.join(HP_dir, 'params_' + str(trial.number) + '.tsv'), 'w')
    for key, value in params.items():
        HP_out.write(str(key)+"\t" + str(value) + "\n")
    HP_out.close()

    #model = build_model(params, train, input_names=["promoter", "terminator"])
    encoder = build_nemo_encoder()
    D = Dense(750, kernel_initializer=ini)(encoder)
    D = BatchNormalization(momentum=Bnorm_momentum)(D)
    D = GELU(approximate=False)(D)

    D = Dense(dense_dim, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=Bnorm_momentum)(D)
    D = GELU(approximate=False)(D)
    D = Dropout(dropout)(D)
    D = Dense(5)(D) # multi-regression output
    model = Model(inputs = [Input(shape=(6200, 4), name="promoter"), Input(shape=(6200, 4), name="terminator")], outputs = D)

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


    # custiom callback for reporting learning rate and momentum
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
    score = min(result.history['val_loss'])

    return score

study.optimize(objective,
               n_trials = n_trials_per_run, # maximum number of trials
               gc_after_trial = True, # garbage collection after each trial
               catch=(ValueError,  # raised when hyperparameters are incompatible
                      MemoryError) # raised when trial uses too much memory
)








