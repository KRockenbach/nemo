#!/usr/bin/python3

##################################################################################
#
# MIT License
#
# Copyright (c) 2026 Kevin Rockenbach, Agnieszka Golicz
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
title: nemo.py
description: From a corresponding annotation (gff3) and assembly (fasta), setup model inputs and make predictions.
author: Kevin Rockenbach
email: kevin.rockenbach@ag.uni-giessen.de
date: 2026-08-07
version: 1.0.0
usage: python nemo.py [-h] -f FASTA -g GFF [-p PREFIX] [-t TMPDIR] [-b BATCH_SIZE] [--weights_ WEIGHTS_]
notes: external dependencies: AGAT, samtools, bedtools
       python dependencies: numpy, pandas, pyarrow, Bio, tensorflow v2.21
       for GPU acceleration: nvidia driver >= 525; CUDA >= 12.5
=========================================================================================================
'''

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from warnings import filterwarnings
filterwarnings("ignore")
from argparse import ArgumentParser
from random import random, seed
from tempfile import TemporaryDirectory
from numpy import column_stack, zeros, round, clip
from pandas import read_csv, read_table, read_feather, Series, merge, DataFrame
from pandas.errors import EmptyDataError
from pickle import load, dump
from Bio import SeqIO
from tensorflow.keras.models import load_model
from math import ceil, floor
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Conv1D, AveragePooling1D, Flatten, Concatenate, Dense, BatchNormalization, Dropout, Activation


seed(1234)

root = os.path.abspath(os.path.dirname(__file__))
weight_path = os.path.join(root, "nemo_v1.0_weights.h5")
scaler_mean = 0.3563692944443251
scaler_scale = 1.0498368135592189

def get_promoter_terminator(outside, inside, index, annotation, outdir):
    prom_out=os.path.join(outdir, "promoters.bed")
    term_out=os.path.join(outdir, "terminators.bed")

    prom_up_part_out=os.path.join(outdir, "promoter_upstream_partials.lst")
    prom_down_part_out=os.path.join(outdir, "promoter_downstream_partials.lst")
    term_up_part_out=os.path.join(outdir, "terminator_upstream_partials.lst")
    term_down_part_out=os.path.join(outdir, "terminator_downstream_partials.lst")

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
        if(l_arr[2] == "mRNA"):
            ch_name = str(l_arr[0])
            mst=l_arr[6] #  mst = mRNA-strand
            ms=int(l_arr[3]) # ms = mRNA-start
            me=int(l_arr[4]) # me = mRNA-end
            attributes=l_arr[8].split(";")
            feature_name={item.split("=")[0]: item.split("=")[1] for item in attributes}["ID"]
            try:
                if (ms > me) or (ms <= 0) or (me > ch_d[ch_name]):
                    print(f"inconsistent bounds for {feature_name}; skipping.")
                    continue
            except KeyError:
                print(f"WARNING: Inconsistent naming of contigs between annotation and assembly. Skipping {feature_name} on {ch_name}.")
                continue

            ############### FORWARD STRAND ################
            if(mst=="+"):
                # CHECK PROMOTER
                start=(ms-outside-1)
                end=(ms+inside-1)
                if start < 0: # upstream sequence of promoter not fully contained
                    prom_up_part.write(feature_name+"\t"+str(abs(start))+"\n")
                    start=0
                if end > ch_d[ch_name]: # downstream sequence of promoter not fully contained
                    prom_down_part.write(feature_name+"\t"+str(end-ch_d[ch_name])+"\n")
                    end=ch_d[ch_name]
                p=[ch_name, str(start), str(end), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")

                # CHECK TERMINATOR
                start=(me-inside)
                end=(me+outside)
                if end > ch_d[ch_name]: # downstream sequence of terminator not fully contained
                    term_down_part.write(feature_name+"\t"+str(end-ch_d[ch_name])+"\n")
                    end=ch_d[ch_name]
                if start < 0: # upstream sequence of terminator not fully contained
                    term_up_part.write(feature_name+"\t"+str(abs(start))+"\n")
                    start=0
                t=[ch_name, str(start), str(end), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")

            ################ REVERSE STRAND ##################
            elif(mst=="-"):
                # CHECK PROMOTER
                start=(me-inside)
                end=(me+outside)
                if end > ch_d[ch_name]: # upstream sequence of promoter not fully contained
                    prom_up_part.write(feature_name+"\t"+str(end-ch_d[ch_name])+"\n")
                    end=ch_d[ch_name]
                if start < 0: # downstream sequence of promoter not fully contained
                    prom_down_part.write(feature_name+"\t"+str(abs(start))+"\n")
                    start=0
                p=[ch_name, str(start), str(end), feature_name+"-promoter", str(0), mst]
                prom.write("\t".join(p)+"\n")

                # CHECK TERMINATOR
                start=(ms-outside-1)
                end=(ms+inside-1)
                if start < 0: # downstream sequence of terminator not fully contained
                    term_down_part.write(feature_name+"\t"+str(abs(start))+"\n")
                    start=0
                if end > ch_d[ch_name]: # upstream sequence of terminator not fully contained
                    term_up_part.write(feature_name+"\t"+str(end-ch_d[ch_name])+"\n")
                    end=ch_d[ch_name]
                t=[ch_name, str(start), str(end), feature_name+"-terminator", str(0), mst]
                term.write("\t".join(t)+"\n")
                #####################
    prom.close()
    term.close()
    prom_up_part.close()
    prom_down_part.close()
    term_up_part.close()
    term_down_part.close()


def merge_data(datadir, prom_path, term_path):

    h=["GENEID", "PROMOTER", "TERMINATOR"]
    merged_data_path = os.path.join(datadir, 'merged.data')
    f = open(merged_data_path, 'w')
    f.write('\t'.join(h) + '\n')
    f.close()

    #create key files containing name and id
    gene_id_path = os.path.join(datadir, 'gene.id.key')
    ginf=open(gene_id_path, 'w') #gene key
    ginf.close()

    gin=1 # gene index (numpy can only andle numerical data)

    # paths to list containing names of partial sequences (upstream or downstream)
    prom_up_path=os.path.join(datadir, "promoter_upstream_partials.lst")
    prom_down_path=os.path.join(datadir, "promoter_downstream_partials.lst")
    term_up_path=os.path.join(datadir, "terminator_upstream_partials.lst")
    term_down_path=os.path.join(datadir, "terminator_downstream_partials.lst")

    try:
        prom_up_lst=read_csv(prom_up_path,
                             header=None, index_col=0, sep="\t")
    except EmptyDataError:
        prom_up_lst=None

    try:
        prom_down_lst=read_csv(prom_down_path,
                               header=None, index_col=0, sep="\t")
    except EmptyDataError:
        prom_down_lst=None

    try:
        term_up_lst=read_csv(term_up_path,
                             header=None, index_col=0, sep="\t")
    except EmptyDataError:
        term_up_lst=None

    try:
        term_down_lst=read_csv(term_down_path,
                               header=None, index_col=0, sep="\t")
    except EmptyDataError:
        term_down_lst=None


    for prom, term in zip(SeqIO.parse(prom_path, "fasta"), SeqIO.parse(term_path, "fasta")):
        sid=prom.id.split("-promoter")[0] #sequence ID
        if prom_up_lst is not None:
            if sid in prom_up_lst.index:
                # pad promoter with Ns from the left
                prom.seq = int(prom_up_lst.loc[sid].item())*"N" + prom.seq
        if term_up_lst is not None:
            if sid in term_up_lst.index:
                # pad terminator with Ns from the left
                term.seq = int(term_up_lst.loc[sid].item())*"N" + term.seq
        if term_down_lst is not None:
            if sid in term_down_lst.index:
                # pad terminator with Ns from the right
                term.seq = term.seq + int(term_down_lst.loc[sid].item())*"N"
        if prom_down_lst is not None:
            if sid in prom_down_lst.index:
                # pad promoter with Ns from the right
                prom.seq = prom.seq + int(prom_down_lst.loc[sid].item())*"N"

        if len(prom.seq) != 6200 or len(term.seq) != 6200:
            print(f"WARNING: Unable to process {sid}. Check inputs. skipping")
            continue

        p = [str(prom.seq), str(term.seq)]
        f = open(merged_data_path, 'a')
        f.write(str(gin)+'\t'+ '\t'.join(p) + '\n')
        f.close()
        ginf=open(gene_id_path, 'a') #gene key
        ginf.write(sid+"\t" + str(gin) + "\n")
        ginf.close()
        gin=gin+1



def get_data(data_file):
    table = read_table(data_file, index_col=0) # read data, gene ID used as index.
    table = table.loc[:,["PROMOTER","TERMINATOR"]]
    assert (not table.isnull().any().any())
    return table

def setup(data_file, outdir):
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    table = get_data(data_file)
    table.to_feather(os.path.join(outdir, f"full.feather"))


def prep(outdir, fasta, gff):
    os.makedirs(outdir, exist_ok = True)
    fasta_path = os.path.join(outdir, "genome.fasta")
    gff_path = os.path.join(outdir, "genes.gff")
    os.system(f"cp {fasta} {fasta_path}")
    os.system(f"cp {gff} {gff_path}")
    rnr=random()
    tmp_gff = os.path.join(outdir, ("tmp_" + str(rnr) + ".gff"))
    # filter and sanitize annotation
    os.system(
        f"""awk 'toupper($3) ~/(MRNA|CDS)/ {{print $0}}' {gff_path} > {tmp_gff} && \
        cat {tmp_gff} > {gff_path} && \
        rm {tmp_gff} || \
        echo "Filtering of GFF failed!" """
    )

    # convert to single line fasta
    tmp_fasta = os.path.join(outdir, ("tmp_" + str(rnr) + ".fa"))
    os.system(
        f"""awk '/^>/ {{printf("\\n%s\\n",$0);next}} {{printf("%s",$0)}}  END {{printf("\\n")}}' < {fasta_path} > {tmp_fasta} && \
        tail -n +2 {tmp_fasta} > {fasta_path} && \
        rm {tmp_fasta} || \
        echo "Single-line conversion failed!" """
    )
    fai_path= (fasta_path + ".fai")
    os.system(f"samtools faidx -n 0 -o {fai_path} {fasta_path}")
    cds_gff_path=gff_path.replace(".gff", ".CDS.gff")
    os.system(f"""awk '$3=="CDS" {{print $0}}' {gff_path} >{cds_gff_path}""")
    print("Masking annotated CDS sequences")
    os.system(f"touch {fai_path}")
    os.system(f"""bedtools maskfasta -fi {fasta_path} -bed {cds_gff_path} -fo {tmp_fasta} 1>{os.path.join(outdir, "bedtools.maskfasta.std")} 2>{os.path.join(outdir, "bedtools.maskfasta.err")} && \
        cat {tmp_fasta} > {fasta_path} && \
        rm {tmp_fasta} || \
        echo "CDS masking failed!" """)

    print("Extracting sequence intervals")
    prom_bed=os.path.join(outdir, "promoters.bed")
    term_bed=os.path.join(outdir, "terminators.bed")
    get_promoter_terminator(5000, 1200, fai_path, gff_path, outdir)
    prom=prom_bed.replace(".bed", ".fasta")
    term=term_bed.replace(".bed", ".fasta")
    os.system(f"touch {fai_path}")
    os.system(f"bedtools getfasta -fi {fasta_path} -bed {prom_bed} -name -s | sed 's/(-)//' | sed 's/(+)//' 1>{prom} 2>{os.path.join(outdir,'bedtools.getfasta.term.err')}")
    os.system(f"bedtools getfasta -fi {fasta_path} -bed {term_bed} -name -s | sed 's/(-)//' | sed 's/(+)//' 1>{term} 2>{os.path.join(outdir,'bedtools.getfasta.prom.err')}")
    merge_data(outdir, prom, term)
    merged_data = os.path.join(outdir, "merged.data")
    os.system(f'gzip {merged_data}')
    print("Setting up model inputs")
    setup(f"{merged_data}.gz", outdir)

def build_nemo():
    promoter = Input(shape=(6200, 4), name="promoter")
    terminator = Input(shape=(6200, 4), name="terminator")

    ini = "glorot_normal"
    bnm = 0.81669
    def get_stride(fraction, size):
        return max([1, ceil(fraction*size)])

    # promoter branch
    P = Conv1D(288, 5, padding = 'same', kernel_initializer = ini)(promoter)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)

    P = Conv1D(106, 8, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)
    P = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(P)

    P = Conv1D(227, 9, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)

    P = Conv1D(187, 10, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)

    P = Conv1D(64, 46, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)
    P = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(P)

    P = Conv1D(152, 62, padding = 'same', kernel_initializer = ini)(P)
    P = BatchNormalization(momentum=bnm)(P)
    P = Activation('gelu')(P)
    P = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(P)

    P = Flatten()(P)

    # terminator branch
    T = Conv1D(219, 4, padding = 'same', kernel_initializer = ini)(terminator)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)

    T = Conv1D(449, 4, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)
    T = AveragePooling1D(21, strides = get_stride(0.7, 21), padding="same")(T)

    T = Conv1D(259, 16, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)

    T = Conv1D(154, 15, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)

    T = Conv1D(105, 20, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)
    T = AveragePooling1D(18, strides = get_stride(0.7, 18), padding="same")(T)

    T = Conv1D(152, 22, padding = 'same', kernel_initializer = ini)(T)
    T = BatchNormalization(momentum=bnm)(T)
    T = Activation('gelu')(T)
    T = AveragePooling1D(8, strides = get_stride(0.7, 8), padding="same")(T)

    T = Flatten()(T)

    # dense layers
    D = Concatenate(axis=1)([P,T])
    D = Dense(750, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=bnm)(D)
    D = Activation('gelu')(D)

    D = Dense(3, kernel_initializer=ini)(D)
    D = BatchNormalization(momentum=bnm)(D)
    D = Activation('gelu')(D)
    D = Dropout(0.00124)(D)
    D = Dense(1)(D) # single-regression output
    return Model(inputs = [promoter, terminator], outputs = D)


def one_hot(seq):
    """ Takes in set of fasta sequences and converts them into
    3D array of one-hot encodings"""
    seq_len = len(seq.iloc[0])
    seqindex = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}
    seq_vec = zeros((seq.count(),seq_len,4), dtype='bool')
    for i in range(seq.count()):
        thisseq = seq.iloc[i]
        for j in range(seq_len):
            try:
                seq_vec[i,j,seqindex[thisseq[j]]] = 1
            except:
                # N = zero column
                pass
    return seq_vec


def get_set4pred(datadir):
    data_dict = {}
    path = os.path.join(datadir, 'full.feather')
    table=read_feather(path)
    data_dict["promoter"] = one_hot(table.loc[:,'PROMOTER'].str.slice(start=0,stop=6200))
    data_dict["terminator"] = one_hot(table.loc[:,'TERMINATOR'].str.slice(start=0,stop=6200))
    data_dict["ID"] = Series(table.index).to_numpy(dtype=int)
    return(data_dict)


def inverse_transform(z):
    n_outputs = 1
    z = (z*scaler_scale) + scaler_mean
    tpm = ((10**z) - 0.1)
    tpm = round(tpm, decimals=1)
    tpm = clip(tpm, a_min=0, a_max=None)
    return tpm


def translate_IDs(ID, datadir, verbose=False):
    # ID is a numpy array with dimensions (#samples, None)
    ID_series = Series(ID[:,None].flatten())
    if verbose:
        print(f"ID shape: {str(ID_series.shape)}")
    ID_series.name = "ID"
    key_file = os.path.join(datadir, "gene.id.key")
    key_df = read_table(key_file, index_col=False, header=None)
    if verbose:
        print(f"key_df shape: {str(key_df.shape)}")
    key_df.columns = ["gene_name", "ID"]
    for ID in ID_series:
        if ID not in key_df.loc[:,"ID"]:
            print(ID)
    # do inner join, preserving the order of the IDs
    merged_df = merge(ID_series, key_df, how="inner", on="ID")
    if verbose:
        print(f"merged_df shape: {str(merged_df.shape)}")
    gene_names = merged_df.loc[:,"gene_name"].to_list()
    return gene_names


def predict(datadir, model_file, prefix, batch):
    modelname = "nemo"
    outdir = datadir

    # load test data
    test = get_set4pred(datadir=datadir)
    inputs = [test["promoter"], test["terminator"]]
    ID = test["ID"]
    print("Loading model")
    model = build_nemo()
    model.load_weights(model_file)

    # get regression predictions
    print("Making predictions")
    preds = model.predict(inputs, batch_size=batch)
    gene_names = translate_IDs(test["ID"], datadir)
    y = inverse_transform(preds)#, scaler)
    mat = column_stack((gene_names, y))
    colnames = ["ID", "Median_TPM"]
    df = DataFrame(mat, columns=colnames)
    f_out = f'{prefix}predictions.tsv'
    df.to_csv(f_out, index=False, header=True, sep='\t')
    print(f"saved predicted values to {f_out}")


def main():
    '''
    From a corresponding annotation (gff3) and assembly (fasta), setup model inputs and make predictions.
    '''
    parser = ArgumentParser(
                    prog='python nemo.py',
                    description='The napus expression model (nemo) predicts median gene expression across tissues under control conditions in Brassica napus. Predictions are saved under <prefix>predictions.tsv',
                    epilog='')
    parser.add_argument('-f', '--fasta', help="Genome assembly in fasta format", required=True)
    parser.add_argument('-g', '--gff', help="Genome annotation in gff format", required=True)
    parser.add_argument('-p', '--prefix', help="Prefix for output file (optional)", default="./")
    parser.add_argument('-t', '--tmpdir', default=root, help="Directory for temporary files (optional)")
    parser.add_argument('-b', '--batch_size' , default="120", help="Batch size used for making predictions (optional)")
    parser.add_argument('--weights_', default=weight_path, help="Path to model weights (optional). By default expected to be in same directory as nemo.py")

    args = parser.parse_args()

    if args.tmpdir:
        os.makedirs(args.tmpdir, exist_ok=True)
        with TemporaryDirectory(dir=args.tmpdir) as tmpdirname:
            prep(tmpdirname, args.fasta, args.gff)
            predict(tmpdirname, args.weights_, args.prefix, int(args.batch_size))
    else:
        with TemporaryDirectory() as tmpdirname:
            prep(tmpdirname, args.fasta, args.gff)
            predict(tmpdirname, args.weights_, args.prefix, int(args.batch_size))

if __name__ == "__main__":
    main()
