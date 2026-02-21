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
title: attrib_utils.py
description: helper functions for attribution and importance related tasks
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      (1) run a given script (e.g. script.py) as: python -m script
      (2) within script, import functions using: from ..utils.attrib_utils import *
=========================================================================================================
'''


import sys, h5py, os, random
import numpy as np
import shap
from tqdm import tqdm
from deeplift.dinuc_shuffle  import *
from random import sample, choice
if __name__ == "main":
    from model_utils import *
else:
    from ..utils.model_utils import *



####################################################
def split_seq(x):
    # split sequences into subsequences contianing only Ns or only called nucleobases
    frag_lst = []
    for seq_idx in range(x.shape[1]):
        current_col = x[:,seq_idx,:]
        current_col = current_col[np.newaxis, :, :] # make array 3D again
        if seq_idx == 0: # first position in sequence
            fragment = current_col
            colsum = np.sum(current_col) # 0 if N, 1 if base
        elif seq_idx == (x.shape[1] - 1): #last position in sequence
            if np.sum(current_col) == colsum:
                fragment = np.append(fragment, current_col, axis=1)
                # store last fragment
                frag_lst.append(fragment)
            else:
                # phase change on last position
                # store previous fragment
                frag_lst.append(fragment)
                # store last position as a single-base fragment
                frag_lst.append(current_col)
        else: #position within seqeunce
            if np.sum(current_col) == colsum:
                fragment = np.append(fragment, current_col, axis=1)
            else:
                # phase change in sequence
                colsum = np.sum(current_col)
                # store fragment
                frag_lst.append(fragment)
                # initialize new fragment
                fragment = current_col
    return frag_lst


############################
def get_dynamic_reference(example, input_type, num_shufs=10, r=None):
    # random number generator for dinuc_shuf can be set as r = np.random.RandomState(1234)
    # example is a list of inputs for a multi-input model. Can also be a single input.
    # input_type is a list of strings that defines the type of each input, either "seq" or "num".
    # numerical features are scaled to a pseudo standard normal distribution (mean=0, sd=1),
    # therefore reference values are randomly drawn from a standard normal distribution
    if not isinstance(example, list):
        example = [example]
    if not isinstance(input_type, list):
        input_type = [input_type]

    references = []
    for idx, input in enumerate(example):
        if input_type[idx] == "seq":
            fragment_lst = split_seq(input)
            shuff_frag_lst = []
            for frag in fragment_lst:
                if np.sum(frag) == 0 or frag.shape[1] == 1: # fragment containing Ns or single-nucleotide fragment
                    shuff_frag_lst.append(np.tile(frag, (num_shufs, 1, 1))) # repeat array num_shufs times
                    #seq_vec = np.zeros((num_shufs,frag.shape[1],4), dtype='bool')
                    #for i in range(num_shufs): # iterate over background
                    #    seqindex = {'A':0, 'C':1, 'G':2, 'T':3}
                    #    rand_seq = "".join(random.choices(["N","A","T","G","C"], k = frag.shape[1]))
                    #    for j in range(frag.shape[1]): # iterate over positions within sequence
                    #        try:
                    #            seq_vec[i,j,seqindex[rand_seq[j]]] = 1 # set corresponding position to 1
                    #        except KeyError: # for N leave position at zero
                    #            pass
                    #shuff_frag_lst.append(seq_vec)
                else: # one-hot array
                    try:
                        shuff_frag_lst.append(dinuc_shuffle(np.squeeze(frag), num_shufs=num_shufs, rng=r))
                    except:
                        print(frag)
                        exit(1)
            references.append(np.concatenate(tuple(shuff_frag_lst), axis=1))
        elif input_type[idx] == "num":
            references.append(r.normal(loc=0.0, scale=1.0, size=(num_shufs,np.squeeze(example[idx]).shape[0])))
            #references.append(np.zeros(shape=(num_shufs, example[idx].shape[0]), dtype="float"))
    return references


#############################
def gradshap(model, x, input_type, num_shufs=10, r=None, verbose=False):
    if not isinstance(x, list):
        x = [x]

    n_samples = x[0].shape[0]
    for sample_idx in range(n_samples):
        # gradient explainer takes in set of background data.
        # to use dynamically generated background data for each sample,
        # attributions need to be calculated for eacha sample individually
        example = [np.expand_dims(ex_input[sample_idx,...], axis=0) for ex_input in x]
        background = get_dynamic_reference(example, input_type,  num_shufs=num_shufs, r=r)
        #explainer
        grad_explainer = shap.GradientExplainer(model, background)

        example_attr = grad_explainer.shap_values(example, nsamples=num_shufs)
        example_attr = example_attr[0]
        #shap output is list of lists (outer list = outputs; inner list = inputs)
        if sample_idx == 0:
            attribs = example_attr
            if verbose:
                print("first sample processed")
        else:
            # list comprehension over all inputs
            attribs = [np.append(attribs[idx], example_attr[idx], axis=0) for idx in range(len(attribs))]
            if verbose:
                if sample_idx > 0 and sample_idx % 50 == 0:
                    print(str(sample_idx) + " samples processed")
    return attribs


########################
def trim_array(in_arr, attr_arr, up, down, start=(-1000), stop=1000):
    trimmed_in_arr = in_arr[:,(up+start):(up+stop),:]
    trimmed_attr_arr = attr_arr[:,(up+start):(up+stop),:]
    assert trimmed_in_arr.shape == trimmed_attr_arr.shape
    assert trimmed_in_arr.shape[1] == (stop-start)
    return trimmed_in_arr, trimmed_attr_arr

############################
def subset_array(attr_array: np.ndarray, seq_array: np.ndarray, name_path: str, group_path_lst: list, get_rest: bool=True, get_names: bool=False) -> np.ndarray:
    import pandas as pd

    def sub_iteration(attr, seq, names, group_idx):
        # get ingroup subset
        group_attr = np.take(attr, group_idx, 0)
        group_seq = np.take(seq, group_idx, 0)
        group_names = [names[idx] for idx in group_idx]
        # get outgroup subset
        rest_idx = [idx for idx in range(attr.shape[0]) if idx not in group_idx]
        rest_attr = np.take(attr, rest_idx, 0)
        rest_seq = np.take(seq, rest_idx, 0)
        rest_names = [names[idx] for idx in rest_idx]
        return group_names, group_seq, group_attr, rest_seq, rest_attr

    names = []
    rest_names = []
    group_idx = []
    rest_seq_lst = []
    rest_attr_lst = []

    for i, group_path in enumerate(group_path_lst):
        group = pd.read_csv(group_path, delimiter=None, header=None).iloc[:,0].tolist()
        group = [x.split(".")[0] for x in group]
        if i == 0:
            names = pd.read_csv(name_path, delimiter=None, header=None).iloc[:,0].tolist()
            group_idx = [names.index(gene) for gene in names if gene in group]
            names, seq, attr, rest_seq, rest_attr = sub_iteration(attr_array, seq_array, names, group_idx)
            rest_seq_lst.append(rest_seq)
            rest_attr_lst.append(rest_attr)
        else:
            group_idx = [names.index(gene) for gene in names if gene in group]
            names, seq, attr, rest_seq, rest_attr = sub_iteration(attr, seq, names, group_idx)
            rest_seq_lst.append(rest_seq)
            rest_attr_lst.append(rest_attr)

    rest_seq = np.concatenate(rest_seq_lst, axis=0)
    rest_attr = np.concatenate(rest_attr_lst, axis=0)

    if get_rest and get_names:
        return names, seq, attr, rest_seq, rest_attr
    elif get_rest:
        return seq, attr, rest_seq, rest_attr
    elif get_names:
        return names, seq, attr
    else:
        return seq, attr


def actual_mean(attr_array: np.ndarray, seq_array: np.ndarray) -> np.ndarray:
    actual_attr = np.multiply(attr_array, seq_array)
    # at each position for each base, only actual attribtions are considered, that correspond to a base actually present at that position
    # sum of attributions at specific position and base, divided by the number of times that base was present at that position
    return(np.sum(actual_attr, axis=0) / np.sum(seq_array, axis=0))


def get_importance(attr_array: np.ndarray, seq_array: np.ndarray) -> np.ndarray:
    # global importance across samples
    # absolute attributions summed over bases, averaged over (non-masked) samples
    importance = np.absolute(attr_array).sum(axis=2).sum(axis=0) / np.sum(seq_array, axis=2).sum(axis=0)
    assert len(importance.shape)==1
    return importance


def get_importance_mean_sd(attr_array: np.ndarray, seq_array: np.ndarray) -> np.ndarray:
    # mean and standard deviation at each position of the attribution array, whith only non-masked bases being taken into account
    mean = np.zeros(shape=attr_array.shape[1])
    sd = np.zeros(shape=attr_array.shape[1])
    seq_base_sum = seq_array.sum(axis=2)
    for p in range(attr_array.shape[1]):
        idx = [True if value == 1 else False for value in seq_base_sum[:,p].tolist()]
        sub =  np.absolute(attr_array[idx,p,:][:,None,:]).sum(axis=2)
        sd[p] = sub.std(axis=0)
        mean[p] = sub.mean(axis=0)
    assert len(mean.shape)==1
    assert len(sd.shape)==1
    return mean, sd


#############################
def save_attribs(input_names, inputs, attr, outdir, method="GradientExplainer"):
    for idx, name in enumerate(input_names):
        if name in ["promoter", "terminator"]:
            seq_arr = inputs[idx]
            shap_arr = attr[idx]
            seq_filename = "".join([name, "_seqs_", method, ".npz"])
            shap_filename = "".join([name, "_shap_", method, ".npz"])
            np.savez_compressed(os.path.join(outdir, seq_filename), seq_arr)
            np.savez_compressed(os.path.join(outdir, shap_filename), shap_arr)
        else: #name = "halflife"
            val_arr = inputs[idx] # axes are not rearranged
            shap_arr = attr[idx]
            val_filename = "".join([name, "_values_", method, ".npz"])
            shap_filename = "".join([name, "_shap_", method, ".npz"])
            np.savez_compressed(os.path.join(outdir, val_filename), val_arr)
            np.savez_compressed(os.path.join(outdir, shap_filename), shap_arr)


def save_clear_inputs(names, inputs, outdir):
    for idx, name in enumerate(names):
        if name in ["promoter", "terminator"]:
            seq_arr = inputs[idx]
            method = "GradientExplainer"
            seq_filename = "".join([name, "_clear_seqs_", method, ".full.npz"])
            np.savez_compressed(os.path.join(outdir, seq_filename), seq_arr)



def modisco_transform(attrib_dir):
    from scipy.signal import savgol_filter
    import pandas as pd
    n_samples = {}
    for seq_name in ["promoter", "terminator"]:
        attr_path = os.path.join(attrib_dir, seq_name + "_shap_GradientExplainer.full.npz")
        in_path = os.path.join(attrib_dir, seq_name + "_seqs_GradientExplainer.full.npz")
        full_attr_arr = np.load(attr_path)['arr_0']
        full_in_arr = np.load(in_path)['arr_0']
        #subset to interval +- 500 bp around TSS/TTS
        if seq_name == "promoter":
            up=5000
            down=1200
        else:
            up=1200
            down=5000

        for start in range(-500, 500, 100):
            in_arr, attr_arr = trim_array(full_in_arr, full_attr_arr, up, down, start=start, stop=(start+100))
            # transform to length last
            in_arr = np.swapaxes(in_arr, 1, 2)
            attr_arr = np.swapaxes(attr_arr, 1, 2)
            keep = []
            seq_len = in_arr.shape[2]
            # masked sequences are not supported by tf-modisco -> use only sequences that are not masked in respective sequence interval
            for sample in range(in_arr.shape[0]):
                if np.sum(in_arr[sample,:,:]) == seq_len:
                    keep.append(True)
                else:
                    keep.append(False)
            in_arr = in_arr[keep,:,:]
            attr_arr = attr_arr[keep,:,:]
            try:
                seq_list = n_samples["sequence"]
                seq_list.append(seq_name)
                n_samples["sequence"] = seq_list
                range_list = n_samples["range"]
                range_list.append(f"{start}to{start+100}")
                n_samples["range"] = range_list
                samples_list = n_samples["num_samples"]
                samples_list.append(in_arr.shape[0])
                n_samples["num_samples"] = samples_list
            except:
                n_samples["sequence"] = [seq_name]
                n_samples["range"] = [f"{start}to{start+100}"]
                n_samples["num_samples"] = [in_arr.shape[0]]

            # save sequences
            out_in = in_path.replace("full", f"modisco_%{start}to{start+100}%")
            np.savez_compressed(out_in, in_arr)
            del in_arr

            # smooth attributions
            magnitude = np.sqrt(np.sum(np.multiply(attr_arr, attr_arr), axis=1, keepdims=True)) # shape = (# examples, 1, seq length)
            unit_vec_arr = np.divide(attr_arr, magnitude)
            np.nan_to_num(unit_vec_arr, copy=False, nan=0.0)
            del attr_arr
            magnitude = savgol_filter(magnitude, window_length = 21, polyorder = 1, axis = 2) #45
            attr_arr = np.multiply(unit_vec_arr, magnitude)
            del unit_vec_arr, magnitude

            # save attributions
            out_attr = attr_path.replace("full", f"modisco_%{start}to{start+100}%")
            np.savez_compressed(out_attr, attr_arr)
            del attr_arr
        del full_in_arr, full_attr_arr
    print("Number of samples based on interval:")
    print(n_samples)
    sample_df = pd.DataFrame.from_dict(n_samples)
    sample_df.to_csv(os.path.join(attrib_dir, "modisco_sample_numbers.tsv"), sep="\t", index=False, header=True)





