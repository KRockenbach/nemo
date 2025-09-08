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
description: converts GFF to BED
author: Agnieszka Golicz
email: agnieszka.golicz@agrar.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage: python nemo/preprocessing/gff2bed.py <GFF> > <BED>
notes: BED file is only needed for gene names
=========================================================================================================
'''

import sys

for l in open(sys.argv[1], 'r'):
    l_arr=l.rstrip().split("\t")
    print(l_arr[0]+"\t"+str(int(l_arr[3])-1)+"\t"+l_arr[4]+"\t"+l_arr[8].split(";")[1].split("=")[1].split(".")[0]+"\t0\t"+l_arr[6])
