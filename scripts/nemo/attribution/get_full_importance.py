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
title: get_full_importance.py
description: calculates importance scores along the entire sequences for all groups in group directory
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python -m nemo.attribution.get_full_importance <group_directory> <distinguishing_feature_between_groups> <model_name> <species>
notes: run within nemo environment
=========================================================================================================
'''


import os, sys
import pandas as pd
import numpy as np
from ..utils.attrib_utils import trim_array, subset_array, get_importance
from ..utils.model_utils import dict_from_tsv, to_bool

group_dir = sys.argv[1]
clustering = sys.argv[2]
model = sys.argv[3]
organism = sys.argv[4]
root=".."
conf_path = os.path.join(root, "model_configs", model, "config.tsv")
config = dict_from_tsv(conf_path)
hp_path = os.path.join(root, "model_configs", model, "hyperparams.tsv")
params = dict_from_tsv(hp_path)
attrib_dir = os.path.join(root, "results", model, organism, "masked_graphpart", "attribs")
name_path = os.path.join(attrib_dir, "gene_names.lst")
well_pred_path = os.path.join(root, "results", model, organism, "masked_graphpart", "IDs", "prediction", "well-predicted.lst")



def export_df(seq_name: str, up: int, down: int) -> None:

    # up and down are the upstream and downstream intervals with respect to the sequence anchor point
    attr_path =os.path.join(attrib_dir, seq_name + "_shap_GradientExplainer.full.npz")
    in_path =os.path.join(attrib_dir, seq_name + "_seqs_GradientExplainer.full.npz")
    attr_arr = np.load(attr_path)['arr_0']
    print(f"Attribution length: {attr_arr.shape[0]}")
    in_arr = np.load(in_path)['arr_0']
    # get subsets

    importance_dict = {}
    if seq_name == "promoter":
        importance_dict["position"] = np.array(range(-5000, 1200))
    else:
        importance_dict["position"] = np.array(range(-1199,5001))

    files = os.listdir(group_dir)
    for group_file in files:
        print(group_file)
        group_path = os.path.join(group_dir, group_file)
        group_descriptor = group_file.split(".")[0]

        if clustering == "prediction":
            group_path_lst=[group_path]
        else:
            group_path_lst=[group_path, well_pred_path]
        group_in_arr, group_attr_arr = subset_array(attr_array=attr_arr, seq_array=in_arr,
                                                    name_path=name_path,
                                                    group_path_lst=group_path_lst,
                                                    get_rest=False, get_names=False)

        group_importance = get_importance(group_attr_arr, group_in_arr) # 1D-array of length 6200
        importance_dict[group_descriptor] = group_importance
        del group_importance

    importance_df = pd.DataFrame(importance_dict)
    importance_path = os.path.join(attrib_dir, seq_name + "_" + clustering + "_importance.tsv")
    importance_df.to_csv(importance_path, sep='\t', header=True, index=False)


promoter = to_bool(config["use_promoter"])
if promoter:
    up = int(params["outside_P"])
    down = int(params["inside_P"])
    export_df(seq_name="promoter", up=up, down=down)
terminator = to_bool(config["use_terminator"])
if terminator:
    up = int(params["inside_T"])
    down = int(params["outside_T"])
    export_df(seq_name="terminator", up=up, down=down)



