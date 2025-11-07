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
title: get_merge_validation.py
description: merges all data needed for prediction into a single file, also sequences intervals that are not fully contained on chromosome are padded
author: Agnieszka Golicz / Kevin Rockenbach
email: agnieszka.golicz@agrar.uni-giessen.de / kevin.rockenbach@ag.uni-giessen.de
date: 2025-11-01
version: 1.0.0
usage:
      python nemo/preprocessing/merge_data.py <data_directory> <FASTA with promoter sequences> <FASTA with terminator sequences>
=========================================================================================================
'''

import sys, os
from Bio import SeqIO
import pandas as pd
import numpy as np


datadir=sys.argv[1]
# set up dictionaries


### TODO make downstream not depend on MEDIAN_EXPRESSION in data set
h=["GENEID", "MEDIAN_EXPRESSION", "PROMOTER", "TERMINATOR"]


###################################
merged_data_path = os.path.join(datadir, 'merged.data')
f = open(merged_data_path, 'w')
f.write('\t'.join(h) + '\n')
f.close()

###################################
# this file should have no header #
###################################
#create key files containing name and id
gene_id_path = os.path.join(datadir, 'gene.id.key')
ginf=open(gene_id_path, 'w') #gene key
ginf.close()
######################################


gin=1 # gene index (numpy can only andle numerical data)


# paths to list containing names of partial sequences (upstream or downstream)
prom_up_path=os.path.join(datadir, "promoter_upstream_partials.lst")
prom_down_path=os.path.join(datadir, "promoter_downstream_partials.lst")
term_up_path=os.path.join(datadir, "terminator_upstream_partials.lst")
term_down_path=os.path.join(datadir, "terminator_downstream_partials.lst")

try:
    prom_up_lst=pd.read_csv(prom_up_path,
                            header=None, index_col=None, sep="\t").iloc[:,0].to_list()
except pd.errors.EmptyDataError:
    print("--- no upstream partial promoter sequences ---")
    prom_up_lst=None
try:
    prom_down_lst=pd.read_csv(prom_down_path,
                          header=None, index_col=None, sep="\t").iloc[:,0].to_list()
except pd.errors.EmptyDataError:
    print("--- no downstream partial promoter sequences ---")
    prom_down_lst=None
try:
    term_up_lst=pd.read_csv(term_up_path,
                        header=None, index_col=None, sep="\t").iloc[:,0].to_list()
except pd.errors.EmptyDataError:
    print("--- no upstream partial terminator sequences ---")
    term_up_lst=None
try:
    term_down_lst=pd.read_csv(term_down_path,
                          header=None, index_col=None, sep="\t").iloc[:,0].to_list()
except pd.errors.EmptyDataError:
    print("--- no downstream partial terminator sequences ---")
    term_down_lst=None


VLKE=0 # number of value list key errors
EDKE=0 # number of exond density key errors
ZDE=0
NEVE=0

prom_path=sys.argv[2]
term_path=sys.argv[3]

transcript_lvl = sys.argv[4]
if transcript_lvl == "True":
    transcript_lvl = True
else:
    transcript_lvl = False

for prom, term in zip(SeqIO.parse(prom_path, "fasta"), SeqIO.parse(term_path, "fasta")):
    sid=prom.id.split("-")[0] #sequence ID
    if not transcript_lvl:
        sid = sid.split(".")[0]
    if prom_up_lst is not None:
        if sid in prom_up_lst:
            # pad promoter with Ns from the left
            prom.seq = (6200 - len(prom.seq))*"N" + prom.seq
            print(f"Padded the promoter-proximal sequence of {sid} from the left")
            if term_up_lst is not None:
                if sid in term_up_lst:
                    # pad terminator with Ns from the left
                    term.seq = (6200 - len(term.seq))*"N" + term.seq
                    print(f"Padded the terminator-proximal sequence of {sid} from the left")
    if term_down_lst is not None:
        if sid in term_down_lst:
            # pad terminator with Ns from the right
            term.seq = term.seq + (6200 - len(term.seq))*"N"
            print(f"Padded the terminator-proximal sequence of {sid} from the right")
            if prom_down_lst is not None:
                if sid in prom_down_lst:
                    # pad promoter with Ns from the right
                    prom.seq = prom.seq + (6200 - len(prom.seq))*"N"
                    print(f"Padded the promoter-proximal sequence of {sid} from the right")



    if len(prom.seq) != 6200 or len(term.seq) != 6200:
        print(f"Problem with {sid}")
        print(f"Promoter length: {str(len(prom.seq))}")
        print(f"Terminator length: {str(len(term.seq))}")
        raise Exception("PaddingError")

    try:

        p = [str(0), str(prom.seq), str(term.seq)]

        f = open(merged_data_path, 'a')
        f.write(str(gin)+'\t'+ '\t'.join(p) + '\n')
        f.close()
        ginf=open(gene_id_path, 'a') #gene key
        ginf.write(sid+"\t" + str(gin) + "\n")
        ginf.close()
        gin=gin+1

    except KeyError: #due to transcript filtering
        VLKE += 1
        continue

    except ZeroDivisionError: # due to CDS length of 0
        ZDE += 1
        continue

    except NoExpressionValueError:
        NEVE += 1
        print(f"No expression values for gene {sid}")
        p = []
        continue


print(f"Number of value list key errors: {str(VLKE)}")
print(f"Number of exon density key errors: {str(EDKE)}")
print(f"Number of ZeroDivisionErrors: {str(ZDE)}")
print(f"Number of NoExpressionValueErrors: {str(NEVE)}")
print("Data has been merged")

