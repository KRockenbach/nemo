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
title: finetue_architecture.py
description: runs third step in hyperparameter optimization
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.hyperparam_tuning.finetune_architecture <name of study> <number trials before periodic restart> <total number of trials to complete>
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


# directory containing training data
datadir = os.path.join("..", "data", "Bnapus", "masked_graphpart_Bn_fold_data")

# CL INPUT 2
study_name = sys.argv[1]

# definition and creation of output directory
# CL INPUT 1
# outdir with config needs do be created beforehand
outdir = os.path.join("..", "model_configs", study_name)

storage = optuna.storages.RDBStorage(
    url="sqlite:///{}/{}.db".format(outdir, study_name),
    heartbeat_interval=30
)

# define study object
n_epo = 5 # max number of epochs
pruner = optuna.pruners.SuccessiveHalvingPruner(min_resource='auto', reduction_factor=4, min_early_stopping_rate=0, bootstrap_count=0)
sampler = TPESampler(consider_prior=True,
                     multivariate=True,
                     group=True,
                     n_startup_trials=100,
                     constant_liar=True) # for distributed optimization

study = optuna.create_study(sampler = sampler,
            study_name = study_name,
            direction = "minimize", # val_loss will be maximized
            pruner = pruner,
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
        summary.write("\t" + "{}:".format(key) + "\t" + "{}".format(value) + "\n")
    summary.close()


    bestHPs = open(os.path.join(outdir, 'hyperparams.tsv'), 'w')
    params = add_manual_hyperparams(bestTrial)
    for key, value in params.items():
        bestHPs.write("{}".format(key)+"\t" + "{}".format(value) + "\n")
    bestHPs.close()

    print("Run Successful")
    print("STUDY STATISTICS" + "\n" + "\n" +
          "Total number of trials:" + "\t" + str(len(study.trials)) + "\n" +
          "Number of pruned trials:" + "\t" + str(n_pruned_trials) + "\n" +
          "Number of failed trials:" + "\t" + str(n_failed_trials) + "\n" +
          "Number of complete trials:" + "\t" + str(n_complete_trials) + "\n")
    exit(0)



def add_manual_hyperparams(trial):
        ###############################
        # Define manual hyperparameters
        ###############################
        params = trial.params

        # shared hyperparameters
        for conv_layer in range(6):
            params[f"pool_size_T_{conv_layer}of6"] = params[f"pool_size_P_{conv_layer}of6"]
            params[f"pool_frac_T_{conv_layer}of6"] = params[f"pool_frac_P_{conv_layer}of6"]
        params["filter_T_5of6"] = params["filter_P_5of6"]

        # fixed hyperparameters
        for branch in ['P', 'T']:
            for layer in range(6):
                params[f"conv_frac_{branch}_{layer}of6"] = 0.0
                params[f"pool_type_{branch}_{layer}of6"] = "avrg"

        params['activation'] = 'gelu'
        params[f'kernel_ini_gelu'] = 'glorot_normal'
        params['use_drop'] = True
        params['use_last_layer_drop'] = True
        params['use_conv_drop'] = False
        params['use_lstm_drop'] = False
        params['ncl_P'] = 6
        params['ncl_T'] = 6
        params['nlstm_P'] = 0
        params['nlstm_T'] = 0
        params['use_exp'] = False
        params['use_Bnorm'] = True
        params['Bnorm_momentum'] = 1 - params['one_minus_Bnorm_momentum']
        params['outside_P'] = 5000
        params['outside_T'] = 5000
        params['inside_P'] = 1200
        params['inside_T'] = 1200
        if params['ndl'] == 2:
            params['dl_drop_0of2'] = 0.0

        return params



##############################################################################
##############################################################################


#Function to build model with trial-specific hyperparameters
def objective(trial):

    K.clear_session() # refresh memory for each trial

    ##########################################################################################
    # define search space
    ##########################################################################################


    batch = trial.suggest_int('batch_size', 16, 32, log=True) # 32

    steps_per_epoch = ceil(n_training / batch)
    print("Number of steps per epoch = " + str(steps_per_epoch))


    use_drop = True
    use_last_layer_drop = True
    use_conv_drop = False
    use_lstm_drop = False

    ncl_P = 6
    ncl_T = 6
    nlstm_P = 0 # number of lstm layers in promoter branch
    nlstm_T = 0 # number of lstm layers in termnator branch

    use_exp = False # mostly useful for getting good first layer representations but can also act as regularization

    #syncronize pooling (size and stride) between branches --> concatenate feature space and feed into LSTM

    for branch in ['P', 'T']:
        # Promoter branch, layer 1
        trial.suggest_int(f"filter_{branch}_0of6", 64, 512, log=True) # 512
        trial.suggest_int(f"kernel_{branch}_0of6", 2, 8, log=True) # 4

        # Promoter branch, layer 2
        trial.suggest_int(f"filter_{branch}_1of6", 64, 512, log=True) # 256
        trial.suggest_int(f"kernel_{branch}_1of6", 2, 8, log=True) # 4

        # Promoter branch, layer 3
        trial.suggest_int(f"filter_{branch}_2of6", 64, 512, log=True) # 128
        trial.suggest_int(f"kernel_{branch}_2of6", 4, 16, log=True) # 16

        # Promoter branch, layer 4
        trial.suggest_int(f"filter_{branch}_3of6", 64, 512, log=True) # 128
        trial.suggest_int(f"kernel_{branch}_3of6", 4, 16, log=True) # 4

        # Promoter branch, layer 5
        trial.suggest_int(f"filter_{branch}_4of6", 64, 512, log=True) # 512
        trial.suggest_int(f"kernel_{branch}_4of6", 16, 64, log=True) # 32

        # Promoter branch, layer 6
        # number of last layer filters shared between branches --> allows for concatenation
        trial.suggest_int(f"kernel_{branch}_5of6", 16, 64, log=True) # 32

    #########################################

    # cross-branch lstm layer
    use_cb_lstm = trial.suggest_categorical('use_cb_lstm', [True, False])
    if use_cb_lstm:
        trial.suggest_int('cb_lstm_units', 16, 128, log=True)

    # Dense layers
    ndl = trial.suggest_int('ndl', 1, 2) # number of dense layers
    trial.suggest_int(f'dl_dim_0of{ndl}', 8, 1024, log=True) # 8
    if ndl == 2:
        trial.suggest_int('dl_dim_1of2', 2, 8, log=True) # --
        # only use dropout on last layer, 0 dropout after first layer, set manually below
        trial.suggest_float('dl_drop_1of2', 5e-3, 5e-1, log=True) # --
    else:
        trial.suggest_float('dl_drop_0of1', 5e-3, 5e-1, log=True) # 0.06766238081567807


    activation = "gelu"
    kernel_ini = "glorot_normal"

    use_Bnorm = True
    one_minus_Bnorm_momentum = trial.suggest_float('one_minus_Bnorm_momentum', 1e-2, 2e-1, log=True)

    # pool size and stride shared between branches
    # hyperparameters in T branch mirrored manually, below
    for conv_layer in range(6):
        trial.suggest_int(f"pool_size_P_{conv_layer}of6", 1, 64, log=True)
        trial.suggest_float(f"pool_frac_P_{conv_layer}of6", 0.5, 1.0, step=0.1)
        if conv_layer == 5:
            trial.suggest_int("filter_P_5of6", 64, 512, log=True)

    params = add_manual_hyperparams(trial)
    print(params)

    ###################################################################################
    # end of search space definition
    ###################################################################################

    ###############################################################################
    # define model
    ###############################################################################

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
    earlystop = EarlyStopping(monitor='val_loss', patience=2, verbose=1, mode='min')

    # fit model
    result = model.fit(
                 [train[i] for i in input_names],
                 train["output"],
                 batch_size = batch,
                 epochs = n_epo,
                 validation_data = ([valid[i] for i in input_names], valid["output"]),
                 callbacks = [pruning, csvlog, term, negR2, earlystop]
    )

    # value used for sampling search space
    score = min(result.history['val_loss'])

    return score

study.optimize(objective,
               n_trials = n_trials_per_run, # maximum number of trials
               gc_after_trial = True, # garbage collection after each trial
#               catch=(ValueError,  # raised when hyperparameters are incompatible
#                      MemoryError) # raised when trial uses too much memory
)








