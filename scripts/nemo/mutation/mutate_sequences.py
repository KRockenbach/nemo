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
title: mutate_sequences.py
description: creates pwm from pfm and motifs from pwm, inserts motif into sequence and makes predictions on mutated sequences
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-02-09
version: 2.0.0
usage: python -m nemo.mutation.mutate_sequences [-o|--organism <species>] [-a|--associated_IDs <path to file containing gene IDs associated with TF family>] [-f|TF_family <TF family name>] [-m|--modelname <name of model>]
notes: run within nemo environment
=========================================================================================================
'''

import sys, os
import numpy as np
import pandas as pd
from ..utils import model_utils as utils #one_hot, build_nemo, dict_from_tsv
from pickle import load
from tqdm import tqdm
from random import sample, choice, shuffle
from argparse import ArgumentParser
import subprocess


n_outputs = 1
out_idx = 2 #median


def inverse_transform(z, scaler, n_outputs = n_outputs, out_idx = out_idx):
        z=z.reshape(-1,n_outputs)
        return (z*scaler.scale_[out_idx])+scaler.mean_[out_idx]


def insert_motifs(prom, term, names, ppm, motif_idx, motif_name, results):
    print(f"mutating")


    def generate_motifs(ppm: pd.DataFrame=ppm) -> list:
        motif_list = []
        scrambled_list = []
        random_list = []

        motif = ""
        random_motif = ""
        for j in range(len(ppm.columns)):
            choices = ["A", "C", "G", "T"]
            weights = ppm.iloc[:,j].to_list()
            motif += np.random.choice(choices, p=weights)
            random_motif += np.random.choice(choices, p=[0.25, 0.25, 0.25, 0.25])
        nucleotides = list(motif)
        shuffle(nucleotides)
        scrambled_motif=''.join(nucleotides)
        return motif, scrambled_motif, random_motif

    def get_sample_idx(array, position) -> np.ndarray:
        sample_bool = np.nonzero(np.sum(array[:,position,:], axis=1))[0] # tuple
        choices = [x for x in range(sample_bool.shape[0]) if sample_bool[x]]
        sample_idx = choice(choices)
        return sample_idx

    def mutate_seq(sample, position, ref_point, motifs, seq_type) -> list:
        abs_position = ref_point + position
        split_list = np.split(sample, [abs_position], axis=0) # axis 0 is length axis after gene selection
        # To keep TSS/TTS coordinate fixed, if position <= TSS/TTS, trim upstream, else trim downstream
        if (abs_position <= ref_point and seq_type == "promoter") or (abs_position < ref_point and seq_type == "terminator"): # trimm upstream (TSS is part of downstream sequence, TTS part of upstream)
            split_list[0] = split_list[0][motifs.shape[1]:,:]
        else: # trim downstream
            split_list[1] = split_list[1][:(split_list[1].shape[0]-motifs.shape[1]),:]
        mutated_list = []
        for m in range(3):
            to_insert = motifs[m,:,:]
            motif_inserted = np.concatenate((split_list[0], to_insert, split_list[1]), axis=0) # concat
            assert(motif_inserted.shape[0]==6200)
            mutated_list.append(np.expand_dims(motif_inserted, axis=0))
        return mutated_list


    mutated_seq_list = [] #
    prom_list = [] #
    term_list = [] #
    rel_idx_lst = [] #
    inserted_motif_lst = [] #
    name_lst = [] #
    motif_type_lst = [] #
    motif_name_lst = [] #


    TSS = 5000
    TTS = 1200
    # every 10th position in +- 500 bp range, offset by motif_idx
    positions = [((-500 + motif_idx) + (10*x)) for x in range(100)]


    for seq_type in ["promoter", "terminator"]:
        for i in tqdm(range(100)):
            motif_types = ["TFBS", "scrambled", "random"]
            type_num = len(motif_types)
            motif_type_lst.extend(motif_types)
            motifs = generate_motifs(ppm)
            inserted_motif_lst.extend(motifs)
            motifs = utils.one_hot(pd.Series(motifs)) # TFBS, scrambled TFBS, random motif
            motif_name_lst.extend([motif_name for j in range(type_num)])
            rel_idx_lst.extend([positions[i] for j in range(type_num)])
            mutated_seq_list.extend([seq_type for j in range(type_num)])

            # draw random sample that is not masked at insertion position
            sample_idx = get_sample_idx(prom, (TSS+positions[i])) if seq_type == "promoter" else get_sample_idx(term, (TTS+positions[i]))
            prom_sample = prom[sample_idx,:,:]
            term_sample = term[sample_idx,:,:]
            name_lst.extend([names[sample_idx] for j in range(type_num)])

            if seq_type == "promoter":
                prom_list.extend(mutate_seq(prom_sample, positions[i], TSS, motifs, seq_type))
                term_list.extend([np.expand_dims(term_sample, axis=0) for j in range(type_num)])
            else:
                term_list.extend(mutate_seq(term_sample, positions[i], TTS, motifs, seq_type))
                prom_list.extend([np.expand_dims(prom_sample, axis=0) for j in range(type_num)])


    df = pd.DataFrame({'name': name_lst, 'motif_type': motif_type_lst, 'motif_name': motif_name_lst, 'motif': inserted_motif_lst,
                       'mutated_sequence_type': mutated_seq_list, 'relative_idx': rel_idx_lst})

    assert(len(prom_list) == len(term_list) == len(df.index))

    prom = np.concatenate(prom_list, axis=0)
    term = np.concatenate(term_list, axis=0)

    del mutated_seq_list #
    del prom_list
    del term_list
    del rel_idx_lst #
    del inserted_motif_lst #
    del name_lst #
    del motif_type_lst #
    del motif_name_lst #

    return df, prom, term



def make_predictions(organism, modelname, config, fold, prom, term):
    print(f"predicting")
    weights = os.path.join("..", "model_weights", modelname, organism, "masked_graphpart")
    model_weights = os.path.join(weights, f"{modelname}_t_{fold}_median.h5")

    if modelname == "nemo":
        model = utils.build_nemo()
    else:
        raise Exception("model not implemented")
    model.load_weights(model_weights)

    input_names = ["promoter", "terminator"]
    inputs = [prom, term]

    # get regression predictions
    batch_size = 32
    mut_preds = model.predict(inputs, batch_size=batch_size)

    # get scaler for inverse transformation
    scaler_path = f"../data/{organism}/masked_graphpart_fold_data/scalers/scaler_{fold}.pkl"
    scaler = load(open(scaler_path, 'rb'))

    # perform inverse transformation
    mut_preds = inverse_transform(mut_preds, scaler)  #scaler expects 2D-array
    return mut_preds


def add_baseline_preds(df, baseline_pred_path):
    pred_df = pd.read_csv(baseline_pred_path, delimiter="\t", header=0, index_col=False)
    baseline = []
    for name in df["name"].tolist():
        baseline.append(pred_df.loc[pred_df["Gene"] == name, "Median_Expression"].iloc[0])
    df["baseline"] = baseline
    return df



def get_ppm(pwm_path: str):
        pwm = pd.read_csv(pwm_path, delimiter="\t", header=None, comment="#")
        ppm = pwm
        for col in range(len(pwm.columns)):
            ppm.iloc[:,col] = (ppm.iloc[:,col]) / (ppm.iloc[:,col].sum())
        return ppm


def exclude_TF_associated(test_set: dict, names: list, TF_family_ID_path: str) -> dict:
    to_exclude = pd.read_csv(TF_family_ID_path, delimiter=None, header=None).iloc[:,0].tolist()
    keep = [False if name in to_exclude else True for name in names]
    for key, value in test_set.items():
        test_set[key] = value[keep] if value.ndim == 1 else value[keep,:] if value.ndim == 2 else value[keep,:,:]
    names=[names[i] for i in range(len(names)) if keep[i]]
    return test_set, names


def main():
    '''
    For every test fold, the corresponding model (trained on the remaining 90% of data) is used to make predictions on mutated genes from the test fold.
    For a given TF-family, 10 random TF motifs per test fold are loaded as PWMs and transformed into PPMs.
    For each PPM, a motif (dynamically generated based on PPM probabilities) is inserted into every 10th position around the TSS and TTS, using the fold index as an offset.
    At the same position, using the same TF-motif, derivative motifs aure also inserted as controls.
    The first control is a scrambled version of the motif, the second is a completely random motif with the same length.
    Background sequences for the insertion are selected from the subset of test set samples that are not masked at the given position, by randomly choosing one of the non-masked sequences at each position.
    Predictions are then made on the mutated samples (one motif insertion per mutated sample, 10 mutated samples per position (+-500 bp) for a given TF family, excluding controls)
    Precalculated baseline predictions made on the origingal sequences are used as reference
    '''
    parser = ArgumentParser(
                    prog='mutate_sequences.py',
                    description='This program inserts motifs into promoter and terminator sequences and predicts expression on the mutated sequences',
                    epilog='')
    parser.add_argument('-o', '--organism', default='Bnapus')
    parser.add_argument('-a', '--associated_IDs')
    parser.add_argument('-f', '--TF_family')
    parser.add_argument('-m', '--modelname', default='nemo')

    args = parser.parse_args()
    TF_family_ID_path = args.associated_IDs
    results = os.path.join("..", "results", args.modelname, args.organism, "masked_graphpart")
    rootdir=".."
    modeldir = os.path.join(rootdir, "model_configs", args.modelname)
    conf_path = os.path.join(modeldir, "config.tsv")
    config = utils.dict_from_tsv(conf_path)

    # TODO IMPORTANT! changed algorithm to take any gene from respective test fold into account, not just moderately expressed!

    outdir = os.path.join(results, "motif_insertion")
    os.makedirs(outdir, exist_ok=True)
    initial = args.organism[0].upper()

    full_df = pd.DataFrame({'name': [], 'motif_type': [], 'motif_name': [], 'motif': [],
                            'mutated_sequence_type': [], 'relative_idx': [],
                            'mutated_pred': [], 'baseline': []})

    # total number of TFs in family
    cmd_string = f"ls ../data/motifs/{args.TF_family} | wc -w"
    num_TFs = int(str(subprocess.check_output(cmd_string, shell=True)).replace("b'", "").replace("\\n'", ""))

    for fold in range(10):
        # random sample of 10 motifs of respectife TF family
        choices = [x+1 for x in range(num_TFs)]
        TF_sample = sample(choices, 10)
        datadir=f"../data/{args.organism}"
        folddir=os.path.join(datadir, "masked_graphpart_fold_data")
        test = utils.get_set(config, outP=5000, inP=1200, outT=5000, inT=1200, datadir=folddir, set="test",
                      test_fold=fold, valid_fold=None)
        test["output"] = test["output"][:,out_idx]
        names = utils.translate_IDs(test["ID"], datadir)

        test, names = exclude_TF_associated(test, names, TF_family_ID_path)

        for motif_idx in range(10):
            pwm_file = str(subprocess.check_output(f"ls ../data/motifs/{args.TF_family} | head -n {TF_sample[motif_idx]} | tail -n 1", shell=True)).replace("b'", "").replace("\\n'", "")
            pwm_path = os.path.join(f"../data/motifs/{args.TF_family}", pwm_file)
            motif_name = pwm_file.split("/")[-1].replace(".pwm", "")
            ppm = get_ppm(pwm_path)
            df, prom, term = insert_motifs(test["promoter"], test["terminator"], names, ppm, motif_idx, motif_name, results)
            mut_preds = make_predictions(args.organism, args.modelname, config, fold, prom, term)

            df["mutated_pred"] = mut_preds

            baseline_pred_path = os.path.join(results, f"nemo{initial}_preds/predictions.t_{fold}.txt")
            df = add_baseline_preds(df, baseline_pred_path)
            full_df = pd.concat([full_df, df], ignore_index=True)
        del test


    f_out = os.path.join(outdir, f"{args.TF_family}.tsv")
    full_df.to_csv(f_out, index=False, header=True, sep='\t')

if __name__ == "__main__":
    main()
