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
title: get_promoter_terminator.py
description: merges all data needed for training into a single file, also sequences intervals that are not fully contained on chromosome are padded
author: Agnieszka Golicz / Kevin Rockenbach
email: agnieszka.golicz@agrar.uni-giessen.de / kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      depending if merging is done for nemo or Xpresso model, either 5 or 9 positional inputs are needed
      python nemo/preprocessing/merge_data.py
      nemo inputs:
                   <data_directory>
                   <TSV with gene names, median expression (only expressed genes)>
                   <FASTA with promoter sequences>
                   <FASTA with terminator sequences>
                   <masking type>
      xpresso inputs:
                   <data_directory>
                   <TSV with gene names, 5P_UTR gc content, 5P_UTR length>
                   <TSV with gene names, 3P_UTR gc content, 3P_UTR length>
                   <TSV with gene names, CDS gc content, CDS length>
                   <TSV with gene names, total gene length, number of CDS features, total exon length>
                   <TSV with gene names, median expression (only expressed genes)>
                   <FASTA with promoter sequences>
                   <FASTA with terminator sequences>
                   <masking type>
=========================================================================================================
'''

import sys, os
from Bio import SeqIO
import pandas as pd
import numpy as np


print(len(sys.argv))
print(sys.argv)
assert(((len(sys.argv) == 10) or (len(sys.argv) == 6)))

datadir=sys.argv[1]
# set up dictionaries

if len(sys.argv) == 10:
    #5UTR
    u5gc_d={} # gc content
    u5len_d={} # length

    for l in open(sys.argv[2]):
        l_arr=l.rstrip().split("\t")
        # gene name used as key
        # same for all subseqent dictionaries
        u5gc_d[l_arr[0].split(".")[0]]=l_arr[1]
        u5len_d[l_arr[0].split(".")[0]]=l_arr[2]


    #3UTR
    u3gc_d={} # gc content
    u3len_d={} # length

    for l in open(sys.argv[3]):
        l_arr=l.rstrip().split("\t")
        u3gc_d[l_arr[0].split(".")[0]]=l_arr[1]
        u3len_d[l_arr[0].split(".")[0]]=l_arr[2]


    #CDS
    cdsgc_d={} # gc content
    cdslen_d={} # length

    for l in open(sys.argv[4]):
        l_arr=l.rstrip().split("\t")
        cdsgc_d[l_arr[0].split(".")[0]]=l_arr[1]
        cdslen_d[l_arr[0].split(".")[0]]=l_arr[2]


    #mRNA/gene model
    cdsex_d={} # number of cds exons
    intronlen_d={} # length
    for l in open(sys.argv[5]):
        l_arr=l.rstrip().split("\t")
        cdsex_d[l_arr[0].split(".")[0]]=l_arr[2]
        gene_len = l_arr[1]
        exon_len = l_arr[3]
        if int(exon_len) > int(gene_len):
            print(l)
            print(f"Exons longer than Gene for {str(l_arr[0].split(".")[0])}")
            continue # annotation error
        intronlen_d[l_arr[0].split(".")[0]]=(int(gene_len) - int(exon_len))



    ##############################
    out_d={}
    ID_col=""
    for idx, l in enumerate(open(sys.argv[6])):
        l_arr=l.rstrip().split("\t")
        out_d[l_arr[0].split(".")[0]]=[str(x) for x in l_arr[1:]]
        if idx == 0:
            ID_col=l_arr[0]
else:
    out_d={}
    ID_col=""
    for idx, l in enumerate(open(sys.argv[2])):
        l_arr=l.rstrip().split("\t")
        out_d[l_arr[0].split(".")[0]]=[str(x) for x in l_arr[1:]]
        if idx == 0:
            ID_col=l_arr[0]


out_names = [ x.upper() for x in out_d[ID_col]]
h=["GENEID"]
h.extend(out_names)
if len(sys.argv) == 10:
    h.extend(["UTR5LEN",
              "CDSLEN",
              "INTRONLEN",
              "UTR3LEN",
              "ORFEXONDENSITY",
              "UTR5GC",
              "CDSGC",
              "UTR3GC",
              "PROMOTER",
              "TERMINATOR"])
    mask_type=sys.argv[9]
else:
    h.extend(["PROMOTER",
              "TERMINATOR"])
    mask_type=sys.argv[5]


###################################
merged_data_path = os.path.join(datadir, 'merged.data.' + mask_type)
f = open(merged_data_path, 'w')
f.write('\t'.join(h) + '\n')
f.close()

###################################
# this file should have no header #
###################################
#create key files containing name and id
gene_id_path = os.path.join(datadir, 'gene.id.key.' + mask_type)
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

if len(sys.argv) == 10:
    prom_path=sys.argv[7]
    term_path=sys.argv[8]
else:
    prom_path=sys.argv[3]
    term_path=sys.argv[4]

for prom, term in zip(SeqIO.parse(prom_path, "fasta"), SeqIO.parse(term_path, "fasta")):
    sid=prom.id.split("-")[0].split(".")[0] #sequence ID
    # Note: Data is matched on gene level, not transcript level.
    if prom_up_lst is not None:
        if sid in prom_up_lst:
            # pad promoter with Ns from the left
            prom.seq = (20000 - len(prom.seq))*"N" + prom.seq
            print(f"Padded the promoter-proximal sequence of {sid} from the left")
            if term_down_lst is not None:
                if sid in term_up_lst:
                    # pad terminator with Ns from the left
                    term.seq = (20000 - len(term.seq))*"N" + term.seq
                    print(f"Padded the terminator-proximal sequence of {sid} from the left")
    if term_down_lst is not None:
        if sid in term_down_lst:
            # pad terminator with Ns from the right
            term.seq = term.seq + (20000 - len(term.seq))*"N"
            print(f"Padded the terminator-proximal sequence of {sid} from the right")
            if prom_down_lst is not None:
                if sid in prom_down_lst:
                    # pad promoter with Ns from the right
                    prom.seq = prom.seq + (20000 - len(prom.seq))*"N"
                    print(f"Padded the promoter-proximal sequence of {sid} from the right")



    if len(prom.seq) != 20000 or len(term.seq) != 20000:
        print(f"Problem with {sid}")
        print(f"Promoter length: {str(len(prom.seq))}")
        print(f"Terminator length: {str(len(term.seq))}")
        raise Exception("PaddingError")

    try:

        p = out_d[sid]
        if len(p) == 0:
            raise NoExpressionValueError(f"No expression values for gene {sid}")
        if len(sys.argv) == 10:
            p.extend([str(u5len_d[sid]),
                str(cdslen_d[sid]),
                str(intronlen_d[sid]), # total intron length
                str(u3len_d[sid]),
                str(round((float(cdsex_d[sid])*1000)/int(cdslen_d[sid]),3)), # exon density
                str(u5gc_d[sid]),
                str(cdsgc_d[sid]),
                str(u3gc_d[sid]),
                str(prom.seq),
                str(term.seq)])
        else:
            p.extend([str(prom.seq),
                      str(term.seq)])

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

