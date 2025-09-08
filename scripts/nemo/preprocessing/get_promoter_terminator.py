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
description: defines intervals around TSS and TTS for sequence extraction from reference, also checks if intervals are fully contained on chromosome
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2025-08-28
version: 1.0.0
usage:
      python nemo/preprocessing/get_promoter_terminator.py <outside_interval_length> <inside_interval_length> <reference.fasta.fai> <annotation.gff> <feature [transcript|mRNA]> <output_directory>
=========================================================================================================
'''

import sys, os

#integer variables for command line arguments (10000  10000)
outside=int(sys.argv[1])
inside=int(sys.argv[2])
index=str(sys.argv[3]) # fasta index file
annotation=str(sys.argv[4]) # genome annotation file
feature=str(sys.argv[5]) # name of the transcript feature (e.g "transcript" or "mRNA", etc.)
outdir=str(sys.argv[6])

# paths to output bed files
prom_out=os.path.join(outdir, "promoters.bed")
term_out=os.path.join(outdir, "terminators.bed")

# paths to list containing names of partial sequences (upstream or downstream)
# partial sequences arise, when the interval extends beyond the boundaries of the chromosome
# partials are to be padded with Ns later (padding is integrated into merge_data.py)
prom_up_part_out=os.path.join(outdir, "promoter_upstream_partials.lst")
prom_down_part_out=os.path.join(outdir, "promoter_downstream_partials.lst")
term_up_part_out=os.path.join(outdir, "terminator_upstream_partials.lst")
term_down_part_out=os.path.join(outdir, "terminator_downstream_partials.lst")

# open output files for writing out in loop
prom=open(prom_out, 'w')
term=open(term_out, 'w')
prom_up_part=open(prom_up_part_out, 'w')
prom_down_part=open(prom_down_part_out, 'w')
term_up_part=open(term_up_part_out, 'w')
term_down_part=open(term_down_part_out, 'w')



ch_d={} # chromosome dictionary
for l in open(index):
    l_arr=l.rstrip().split("\t")

    ch_name = str(l_arr[0])
    ch_len = int(l_arr[1])
    ch_d[ch_name]=ch_len # length of each chromosome

for l in open(annotation, encoding="utf8", errors='ignore'):
    if(l.startswith("#")): # ignore headers
        continue
    l_arr=l.rstrip().split("\t")
    if(l_arr[2]==feature):
        ch_name = str(l_arr[0])
        mst=l_arr[6] #  mst = mRNA-strand
        ms=int(l_arr[3]) # ms = mRNA-start
        me=int(l_arr[4]) # me = mRNA-end
        feature_name=l_arr[8].split(";")[0].split("=")[1].split(".")[0]


        # check if invervals are fully contained on chromosome (GFF is 1-based)
        ############### FORWARD STRAND ################
        if(mst=="+"):

            ##################
            # CHECK PROMOTER #
            ##################
            # 10000 bp long upstream interval starting at first basepair would be [1,10000]
            # for interval to be fully contained on chromosome, gene would have to have start coordinate 10001
            # ---> outside promoter interval not fully contained when start coordinate - interval length is less than 1
            if (ms-outside) < 1: # upstream sequence of promoter not fully contained
                # write name of partial sequence
                prom_up_part.write(feature_name+"\n")
                # write sequence coordinates into bed (starts at 0)
                p=[ch_name, str(0), str(ms+inside-1), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            # for downstream interval of promoter, sum of coordinate and interval length must be subtracted by 1,
            # to not count interval start coordinate twice
            elif (ms+inside-1) > ch_d[ch_name]: # downstream sequence of promoter not fully contained
                # write name of partial sequence
                prom_down_part.write(feature_name+"\n")
                # write sequence coordinates into bed (ends at chromosome end)
                p=[ch_name, str(ms-outside-1), str(ch_d[ch_name]), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            else: # promoter intervals fully contained
                # -1 is for conversion of start coordinate to bed format
                p=[ch_name, str(ms-outside-1), str(ms+inside-1), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            ##################
            ##################



            ####################
            # CHECK TERMINATOR #
            ####################
            # in 1-based closed coordinate system, end coordinate gives interval length from start until and including that coordinate
            # for downstream terminator interval, outside interval lengths just needs to be added ontop
            if (me+outside) > ch_d[ch_name]: # downstream sequence of terminator not fully contained
                # write name of partial sequence
                term_down_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (ends at chromosome end)
                t=[ch_name, str(me-inside), str(ch_d[ch_name]), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            # for inside interval, feature coordinate is part of interval
            # ---> inside terminator interval not fully contained when end coordinate - interval length is less than 0,
            # since gene ending at GFF-coordinate 10000 could have fully contained interval [1,10000]
            elif (me-inside) < 0: # upstream sequence of terminator not fully contained
                # write name of partial sequence
                term_up_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (starts at 0)
                t=[ch_name, str(0), str(me+outside), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            else: # terminator intervals fully contained
                t=[ch_name, str(me-inside), str(me+outside), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            ##################
            ##################


        ################ REVERSE STRAND ##################
        elif(mst=="-"):


            ##################
            # CHECK PROMOTER #
            ##################
            if (me+outside) > ch_d[ch_name]: # upstream sequence of promoter not fully contained
                # write name of partial sequence
                prom_up_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (ends at chromosome end)
                p=[ch_name, str(me-inside), str(ch_d[ch_name]), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            elif (me-inside) < 0: # downstream sequence of promoter not fully contained
                # write name of partial sequence
                prom_down_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (starts at 0)
                p=[ch_name, str(0), str(me+outside), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            else: # promoter intervals fully contained
                p=[ch_name, str(me-inside), str(me+outside), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")
            #####################
            #####################



            ####################
            # CHECK TERMINATOR #
            ####################
            if (ms-outside) < 1: # downstream sequence of terminator not fully contained
                # write name of partial sequence
                term_down_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (starts at 0)
                t=[ch_name, str(0), str(ms+inside-1), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            elif (ms+inside-1) > ch_d[ch_name]: # upstream sequence of terminator not fully contained
                # write name of partial sequence
                term_up_part.write(feature_name+"\n")
                # write sequence coordiantes into bed (ends at chromosome end)
                t=[ch_name, str(ms-outside-1), str(ch_d[ch_name]), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            else: # terminator intervals fully contained
                t=[ch_name, str(ms-outside-1), str(ms+inside-1), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
            #####################
            #####################


prom.close()
term.close()
prom_up_part.close()
prom_down_part.close()
term_up_part.close()
term_down_part.close()

