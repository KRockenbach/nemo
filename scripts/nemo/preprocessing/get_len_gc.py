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
title: get_len_gc.py
description: calculated lengths and GC contents of CDS or UTR features
author: Agnieszka Golicz
email: agnieszka.golicz@agrar.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      python nemo/preprocessing/get_len_gc.py <CDS.fasta|5P_UTR.fasta|3P_UTR.fasta> <promoters.bed>
notes: BED file is only needed for gene names
=========================================================================================================
'''

import sys
from Bio import SeqIO

# a dictionary for counting each base of each sequence ID
a_d={}
c_d={}
g_d={}
t_d={}
n_d={}

# dictionaries for sequence length and gc content
len_d={}
gc_d={}

for seq in SeqIO.parse(sys.argv[1], "fasta"):
    a = seq.seq.upper().count('A')
    c = seq.seq.upper().count('C')
    g = seq.seq.upper().count('G')
    t = seq.seq.upper().count('T')
    n = seq.seq.upper().count('N')

    key = seq.id.split("-")[0].split(":")[0].split(".")[0] #sequence ID
    # Some genes have multi-exon UTRs
    if key in a_d.keys(): # add count to previous entry
        a_d[key] += a
        c_d[key] += c
        g_d[key] += g
        t_d[key] += t
        n_d[key] += n
    else: # create new entry
        a_d[key] = a
        c_d[key] = c
        g_d[key] = g
        t_d[key] = t
        n_d[key] = n

    # calculate gc content based on total counts for given sequence ID
    gc=(float(c_d[key]) + float(g_d[key]))/(float(a_d[key])+float(c_d[key])+float(g_d[key])+float(t_d[key])+float(n_d[key]))

    # update length and gc dictionaries
    gc_d[key]=str(round(gc, 3))

    if key in len_d.keys():
        len_d[key]=str(int(len_d[key])+len(seq))
    else:
        len_d[key]=str(len(seq))


for l in open(sys.argv[2]):
    l_arr=l.rstrip().split("\t")
    key = l_arr[3].split("-")[0].split(":")[0].split(".")[0]
    if(key in gc_d.keys()):
        print(key+"\t"+gc_d[key]+"\t"+len_d[key])
    else: #make sure to also include transcripts that do not have a UTR
        print(key+"\t0\t0")
