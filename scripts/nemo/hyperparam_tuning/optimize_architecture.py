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
title: optimized_architecture.py
description: runs first step in hyperparameter optimization
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.hyperparam_tuning.optimize_architecture <name of study> <number trials before periodic restart> <total number of trials to complete>
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
n_epo = 10 # max number of epochs
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


# get full-length (default) data for the training and validation set
# input sequences will be subsampled within the objective function
#########################
train = get_set(config, datadir=datadir, set="training", test_fold=test_fold, valid_fold=valid_fold)
valid = get_set(config, datadir=datadir, set="validation", test_fold=test_fold, valid_fold=valid_fold)

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
    if best_params["activation"] == "selu":
        best_params["kernel_ini_selu"] = "lecun_normal"
    else:
        activation = best_params["activation"]
        if activation == "relu" or activation == "leaky_relu" or activation == "prelu":
            best_params["kernel_ini_{activation}"] = "he_normal"
        # else initialization already added to param dict by optuna

    for key, value in bestTrial.params.items():
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

    ##########################################################################################
    # define search space
    ##########################################################################################

    # input sequence intervals
    # search space for the inside sequences is restricted to 5 kb to avoid excessive overlap between sequence inputs

    outside = trial.suggest_int('outside', 1000, 15000, step=1000)
    inside = trial.suggest_int('inside', 100, 1500, step=100)

    # create datasets with subsampled sequences
    # full data is loaded into memory once when first running the script
    # subsets need to be created using .copy() method so that the original set is preserved
    sub_train = {}
    sub_valid = {}
    # locations within respective sequence intervals
    TSS=15000
    TTS=5000
    for key in train.keys():
        if key == "promoter":
            sub_train[key] = train[key].copy()[:,TSS-outside:TSS+inside,:]
            sub_valid[key] = valid[key].copy()[:,TSS-outside:TSS+inside,:]
        elif key == "terminator":
            sub_train[key] = train[key].copy()[:,TTS-inside:TTS+outside,:]
            sub_valid[key] = valid[key].copy()[:,TTS-inside:TTS+outside,:]
        else:
            sub_train[key] = train[key]
            sub_valid[key] = valid[key]



    batch_exp = trial.suggest_int('batch_exponent', 2, 7) #4 to 128
    batch = 2**batch_exp

    steps_per_epoch = ceil(n_training / batch)
    print("Number of steps per epoch = " + str(steps_per_epoch))


    use_drop = trial.suggest_categorical("use_drop", [True, False])
    if use_drop:
        use_last_layer_drop = trial.suggest_categorical("use_last_layer_drop", [True, False])
        if use_last_layer_drop:
            use_conv_drop = False
            use_lstm_drop = False
        else:
            use_conv_drop = trial.suggest_categorical("use_conv_drop", [True, False])
            use_lstm_drop = trial.suggest_categorical("use_lstm_drop", [True, False])


    ncl_P = trial.suggest_int('ncl_P', 1, 6) # number of convolutional layers in promoter branch
    ncl_T = trial.suggest_int('ncl_T', 1, 6) # number of convolutional layers in termnator branch
    nlstm_P = trial.suggest_int('nlstm_P', 0, 2) # number of lstm layers in promoter branch
    nlstm_T = trial.suggest_int('nlstm_T', 0, 2) # number of lstm layers in termnator branch

    use_exp = trial.suggest_categorical("use_exp", [True, False]) # mostly useful for getting good first layer representations but can also act as regularization

    for branch in ["P", "T"]:
        ncl = eval(f"ncl_{branch}")
        for layer in range(ncl):
            flt = f'filter_{branch}_{layer}of{ncl}'
            krn = f'kernel_{branch}_{layer}of{ncl}'
            pl_size = f'pool_size_{branch}_{layer}of{ncl}'
            pl_frac = f'pool_frac_{branch}_{layer}of{ncl}'
            pl_type = f'pool_type_{branch}_{layer}of{ncl}'
            cv_frac = f'conv_frac_{branch}_{layer}of{ncl}'
            trial.suggest_int(("exponent_" + flt), 5, 9) # 32 to 512
            if layer == 0 and use_exp:
                trial.suggest_int(("exponent_" + krn + "_ExpActivation"), 4, 5) # 16 to 32 (kernel should be large enough to fit whole motif)
            else:
                trial.suggest_int(("exponent_" + krn), 2, 5) # 4 to 32
            trial.suggest_int(("exponent_" + pl_size), 0, 6) # 1 to 64
            trial.suggest_float(cv_frac, 0.0, 0.5, step=0.1)
            trial.suggest_categorical(pl_type, ["max","avrg"])
            if trial.params[pl_type] in ["max", "avrg"]: # set pool stride
                # overlapping pooling windows could help assemble motifs in deeper layers
                # discrete search space makes full window stride and stride 1 more likely
                trial.suggest_float(pl_frac, 0.5, 1.0, step=0.1)
                # stride is calculated as max(frac*kernel_size, 1),
                # so frac = 0.0 results in stride = 1
                # attention pooling does not have stride parameter
            if use_drop:
                if use_conv_drop:
                    cv_drop = f'cv_drop_{branch}_{layer}of{ncl}'
                    trial.suggest_float(cv_drop, 1e-6, 1.0, log=True)


        nlstm = eval(f"nlstm_{branch}") # number of LSTM layers for individual branch
        if nlstm > 0:
            skip_type = trial.suggest_categorical('skip_type', ['none', 'add', 'multiply', 'concatenate']) # implementation of skip connection
            for layer in range(nlstm):
                lstm_units = f'lstm_units_{branch}_{layer}of{nlstm}'
                if (layer != (nlstm-1)) or (skip_type in ['concatenate', 'none']):
                    trial.suggest_int(("exponent_" + lstm_units), 5, 9) # 32 to 512
                if use_drop:
                    if use_lstm_drop:
                        lstm_drop = f'lstm_drop_{branch}_{layer}of{nlstm}'
                        trial.suggest_float(lstm_drop, 1e-6, 1.0, log=True)


    ndl = trial.suggest_int('ndl', 1, 3) # number of regression dense layers
    for layer in range(ndl):
        dl_dim = f'dl_dim_{layer}of{ndl}'
        trial.suggest_int(("exponent_" + dl_dim), 1, 10)
        if use_drop:
            if not (use_last_layer_drop and layer < (ndl-1)): # else set to zero manually (below)
                dl_drop = f'dl_drop_{layer}of{ndl}'
                trial.suggest_float(dl_drop, 1e-6, 1.0, log=True)



    activation = trial.suggest_categorical('activation', ["relu", "leaky_relu", "selu", "prelu",
                                                          "gelu", "elu", "swish", "softplus", "mish"])
    if activation == "selu":
        kernel_ini = 'lecun_normal'
    else:
        if activation == "relu" or activation == "leaky_relu" or activation == "prelu":
            # relu-like activation functions should have He initialization
            kernel_ini = 'he_normal'
            if activation == "leaky_relu" or activation == "prelu":
                # in case of prelu, alpha is used for initialization
                relu_alpha = trial.suggest_float('relu_alpha', 5e-6, 0.5, log=True)
        else: # gelu, elu, swish, mish, softplus
            # appropriate initialization not as straight forward
            trial.suggest_categorical(f'kernel_ini_{activation}', ['he_normal', 'glorot_normal'])
            if activation == 'elu':
                elu_alpha = trial.suggest_float('elu_alpha', 0.2, 2.0, log=True) # log-centered around 1.0, which is default value


#    use_Bnorm = trial.suggest_categorical("use_Bnorm", [True, False]) # TF2.x implementation not supported by DeepSHAP, but should be supportet by SHAP gradient explainer
    use_Bnorm = True
    if use_Bnorm:
        one_minus_Bnorm_momentum = trial.suggest_float('one_minus_Bnorm_momentum', 1e-4, 4e-1, log=True)



###############################
# Define manual hyperparameters
###############################

    params = trial.params
    # add manually defined hyperparameters to dictionary
    try:
        kernel_ini = params[f'kernel_ini_{activation}']
    except KeyError:
        params[f'kernel_ini_{activation}'] = kernel_ini

    params['batch_size'] = batch

    params['use_Bnorm'] = use_Bnorm
    if use_Bnorm:
        params['Bnorm_momentum'] = 1 - one_minus_Bnorm_momentum

    params['outside_P'] = outside
    params['outside_T'] = outside
    params['inside_P'] = inside
    params['inside_T'] = inside

    if use_drop:
        if use_last_layer_drop:
            params['use_conv_drop'] = False
            params['use_lstm_drop'] = False


    for branch in ["P", "T"]:
        ncl = eval(f"ncl_{branch}")
        for layer in range(ncl):
            flt = f'filter_{branch}_{layer}of{ncl}'
            krn = f'kernel_{branch}_{layer}of{ncl}'
            pl_size = f'pool_size_{branch}_{layer}of{ncl}'
            params[flt] = 2**(int(params[("exponent_" + flt)]))
            if layer == 0 and use_exp:
                params[krn] = 2**(int(params[("exponent_" + krn + "_ExpActivation")]))
            else:
                params[krn] = 2**(int(params[("exponent_" + krn)]))
            params[pl_size] = 2**(int(params[("exponent_" + pl_size)]))
        nlstm = eval(f"nlstm_{branch}")
        for layer in range(nlstm):
            lstm_units = f'lstm_units_{branch}_{layer}of{nlstm}'
            if layer == (nlstm-1) and skip_type in ['add', 'multiply']:
                # dimensions must match for element-wise addition/multiplication
                # bidirectional layer will produce twice the number of features -> duplicate of feature map will be concatenated in the skip connection
                params[lstm_units] = params[f'filter_{branch}_{ncl-1}of{ncl}']
            else:
                params[lstm_units] = 2**(int(params[("exponent_" + lstm_units)]))

    for layer in range(ndl):
        dl_dim = f'dl_dim_{layer}of{ndl}'
        params[dl_dim] = 2**(int(params[("exponent_" + dl_dim)]))
        if use_drop:
            if use_last_layer_drop and layer < (ndl-1): # zero dropout on intermediate layers (to avoid interference with batchnorm)
                dl_drop = f'dl_drop_{layer}of{ndl}'
                params[dl_drop] = 0.0

    '''
    hyperparameter values can be called either through variable name, if assigned,
    or through params['key'] = value
    '''


    ###################################################################################
    # end of search space definition
    ###################################################################################

    ###############################################################################
    # define model
    ###############################################################################

    print(params)

    model = build_model(params, sub_train, input_names=["promoter", "terminator"])

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

    # fit model
    result = model.fit(
                 [sub_train[i] for i in input_names],
                 sub_train["output"],
                 batch_size = batch,
                 epochs = n_epo,
                 validation_data = ([sub_valid[i] for i in input_names], sub_valid["output"]),
                 callbacks = [pruning, csvlog, term, negR2]
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








