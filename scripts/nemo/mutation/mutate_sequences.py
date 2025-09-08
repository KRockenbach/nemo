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
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.mutation.mutate_sequences <path to pwm of TF> <species> <path to file containing gene IDs associated with TF family> <TF family name>
notes: run within nemo environment
=========================================================================================================
'''

import sys, os
import numpy as np
import pandas as pd
from ..utils.model_utils import one_hot, build_nemo
from pickle import load
from tqdm import tqdm
from random import sample, choice


pwm_path = sys.argv[1]
organism = sys.argv[2]
TF_family_ID_path = sys.argv[3]
family = sys.argv[4]
n_insert=20
results = os.path.join("..", "results", "nemo", organism, "masked_graphpart")

high_path = os.path.join(results, "IDs/expression/high_expr_ids.lst")
medium_path = os.path.join(results, "IDs/expression/medium_expr_ids.lst")
low_path = os.path.join(results, "IDs/expression/low_expr_ids.lst")

outdir = os.path.join(results, f"{family}_insertion")
os.makedirs(outdir, exist_ok=True)
initial = organism[0].upper()
original_pred_path = os.path.join(results, f"nemo{initial}_preds/predictions.full.txt")



def insert_motifs(pwm_path, organism, group_path_list, n_insert):

    prom_path = os.path.join(results, "attribs/promoter_seqs_GradientExplainer.full.npz")
    term_path = os.path.join(results, "attribs/terminator_seqs_GradientExplainer.full.npz")
    name_path = os.path.join(results, "attribs/gene_names.lst")

    def subset_array(prom_array, term_array, name_path, group_path_list, TSS, TTS):

        def sub_iteration(prom, term, names, group_idx):
            # get ingroup subset
            group_prom = np.take(prom, group_idx, 0)
            group_term = np.take(term, group_idx, 0)
            group_names = [names[idx] for idx in group_idx]
            return group_names, group_prom, group_term

        names = pd.read_csv(name_path, delimiter=None, header=None).iloc[:,0].tolist()
        prom = prom_array
        term = term_array
        group_idx = []

        for i, group_path in enumerate(group_path_list):
            group = pd.read_csv(group_path, delimiter=None, header=None).iloc[:,0].tolist()
            group = [x.split(".")[0] for x in group]
            group_idx = [names.index(gene) for gene in names if gene not in group] # group gets excluded!
            names, prom, term = sub_iteration(prom, term, names, group_idx)

        prom_pick = np.asarray(prom[:,(TSS-1000):(TSS+1000),:].sum(axis=2).sum(axis=1) >= 20).nonzero()[0].tolist()
        term_pick = np.asarray(term[:,(TTS-1000):(TTS+1000),:].sum(axis=2).sum(axis=1) >= 20).nonzero()[0].tolist()
        pick = [item for item in prom_pick if item in term_pick] # only pick from genes that have at least 20 unmasked bases around TSS and TTS
        keep = sample(pick, 1000) # doing it for all genes at once leads to OOM
        names = [names[i] for i in keep]
        prom = prom[keep,:,:]
        term = term[keep,:,:]
        return names, prom, term


    def generate_motif_list(pwm_path: str=pwm_path, n_motifs: int=n_insert) -> list:
        motif_list = []
        pwm = pd.read_csv(pwm_path, delimiter="\t", header=None, comment="#")
        ppm = pwm
        for col in range(len(pwm.columns)):
            ppm.iloc[:,col] = (ppm.iloc[:,col]) / (ppm.iloc[:,col].sum())
        for i in range(n_motifs):
            motif = ""
            for j in range(len(ppm.columns)):
                choices = ["A", "C", "G", "T"]
                weights = ppm.iloc[:,j].to_list()
                motif += np.random.choice(choices, p=weights)
            motif_list.append(motif)
        return motif_list


    TSS=5000
    TTS=1200
    prom_arr = np.load(prom_path)['arr_0']
    term_arr = np.load(term_path)['arr_0']
    names, prom_arr, term_arr = subset_array(prom_array=prom_arr, term_array=term_arr, name_path=name_path, group_path_list=group_path_list, TSS=TSS, TTS=TTS)
    mutated_prom_list = []
    original_prom_list = []
    mutated_term_list = []
    original_term_list = []

    prom_idx_lst = []
    term_idx_lst = []
    inserted_motif_lst = []
    new_name_lst = []
    parent_name_lst = []
    for i in tqdm(range(len(names))):
        motif_list = generate_motif_list()
        motifs = one_hot(pd.Series(motif_list)) # new set of probabilistically generated motifs for each sequence
        for j in range(n_insert):
            to_insert=motifs[j,:,:]
            parent_name_lst.append(names[i])
            new_name_lst.append(names[i]+f":{j}")
            inserted_motif_lst.append(motif_list[j])
            idx_choices = np.asarray(prom_arr[i,(TSS-1000):(TSS+1000),:].sum(axis=1) == 1).nonzero()[0].tolist()
            random_prom_idx = choice(idx_choices) + TSS - 1000 # samples from -1000 to 999 around TSS, except masked regions
            prom_idx_lst.append(random_prom_idx)
            split_prom_list = np.split(prom_arr[i,:,:], [random_prom_idx], axis=0) # axis 0 is length axis after gene selection
            # To keep TSS coordinate fixed, if random_prom_idx <= TSS, trim upstream, else trim downstream (TSS is part of 5kb downstream sequence!)
            if random_prom_idx <= TSS: # trimm upstream
                split_prom_list[0] = split_prom_list[0][motifs.shape[1]:,:]
            else: # trim downstream
                split_prom_list[1] = split_prom_list[1][:(split_prom_list[1].shape[0]-motifs.shape[1]),:]
            motif_inserted_prom = np.concatenate((split_prom_list[0], to_insert, split_prom_list[1]), axis=0) # concat
            assert(motif_inserted_prom.shape[0]==6200)
            mutated_prom_list.append(np.expand_dims(motif_inserted_prom, axis=0))
            original_prom_list.append(np.expand_dims(prom_arr[i,:,:], axis=0))

            idx_choices = np.asarray(term_arr[i,(TTS-1000):(TTS+1000),:].sum(axis=1) == 1).nonzero()[0].tolist()
            random_term_idx = choice(idx_choices) + TTS - 1000 # samples from -1000 to 999 around TTS, except masked region
            term_idx_lst.append(random_term_idx)
            split_term_list = np.split(term_arr[i,:,:], [random_term_idx], axis=0)
            # To keep TTS coordinate fixed, if random_term_idx < TTS, trim upstream, else trim downstream (TTS is part of 1.2 kb upstream sequence!)
            if random_term_idx < TTS: # trim upstream
                split_term_list[0] = split_term_list[0][motifs.shape[1]:,:]
            else: # trim downstream
                split_term_list[1] = split_term_list[1][:(split_term_list[1].shape[0]-motifs.shape[1]),:]
            motif_inserted_term = np.concatenate((split_term_list[0], to_insert, split_term_list[1]), axis=0)
            assert(motif_inserted_term.shape[0]==6200)
            mutated_term_list.append(np.expand_dims(motif_inserted_term, axis=0))
            original_term_list.append(np.expand_dims(term_arr[i,:,:],axis=0))
    df = pd.DataFrame({'new_name': new_name_lst, 'parent_name': parent_name_lst,
                       'motif': inserted_motif_lst, 'prom_insert_idx': prom_idx_lst, 'term_insert_idx': term_idx_lst})

    assert(len(mutated_prom_list) == len(original_prom_list) == len(mutated_term_list) == len(original_term_list))

    mutated_prom = np.concatenate(mutated_prom_list, axis=0)
    original_prom = np.concatenate(original_prom_list, axis=0)
    mutated_term = np.concatenate(mutated_term_list, axis=0)
    original_term = np.concatenate(original_term_list, axis=0)

    assert mutated_prom.shape[0] == (n_insert*prom_arr.shape[0])
    assert mutated_term.shape[0] == (n_insert*term_arr.shape[0])
    return df, mutated_prom, original_term, original_prom, mutated_term



def make_predictions(organism, mutated_prom, original_prom, mutated_term, original_term):

    weights = os.path.join("..", "model_weights", "nemo", organism, "masked_graphpart")
    model_weights = os.path.join(weights, "nemo_full.h5")

    model = build_nemo()
    model.load_weights(model_weights)
    # get list of input names
    input_names = ["promoter", "terminator"]

    # get list of input data
    prom_mut_inputs = [mutated_prom, original_term]
    term_mut_inputs = [original_prom, mutated_term]

    # get regression predictions
    batch_size=16
    prom_mut_preds = model.predict(prom_mut_inputs, batch_size=16)
    term_mut_preds = model.predict(term_mut_inputs, batch_size=16)

    # get scaler for inverse transformation
    scaler_path = f"../data/{organism}/masked_graphpart_fold_data/scalers/scaler_full.pkl"
    scaler = load(open(scaler_path, 'rb'))

    # perform inverse transformation
    n_outputs=1
    prom_mut_preds = scaler.inverse_transform(prom_mut_preds.reshape(-1,n_outputs))  #scaler expects 2D-array
    term_mut_preds = scaler.inverse_transform(term_mut_preds.reshape(-1,n_outputs))
    return prom_mut_preds, term_mut_preds


def add_original_preds(df, original_pred_path):
    pred_df = pd.read_csv(original_pred_path, delimiter="\t", header=0, index_col=False)
    baseline = []
    for name in df["parent_name"].tolist():
        baseline.append(pred_df.loc[pred_df["Gene"] == name, "Median_Expression"].iloc[0])
    df["baseline"] = baseline
    return df



for exp_group in ["medium"]: #["high", "medium", "low"]:
    if exp_group=="high":
        group_path_list = [medium_path, low_path] # groups to exclude
    elif exp_group=="medium":
        group_path_list = [low_path, high_path]
    else:
        group_path_list = [high_path, medium_path]
    group_path_list.append(TF_family_ID_path)

    df, mutated_prom, original_term, original_prom, mutated_term = insert_motifs(pwm_path, organism, group_path_list, n_insert)
    prom_mut_preds, term_mut_preds = make_predictions(organism, mutated_prom, original_prom, mutated_term, original_term)

    df["prom_mut_preds"] = prom_mut_preds
    df["term_mut_preds"] = term_mut_preds

    df = add_original_preds(df, original_pred_path)

    motif_type = pwm_path.split("/")[-1].replace(".pwm", "")
    f_out = os.path.join(outdir, f"{motif_type}_{exp_group}_exp.tsv")
    df.to_csv(f_out, index=False, header=True, sep='\t')
