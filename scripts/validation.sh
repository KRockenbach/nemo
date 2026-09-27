#!/usr/bin/env bash

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


#========================================================================================================
#title: validation.sh
#description: validates nemo predicitons on real world data
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-11-07
#version: 1.0.0
#usage: bash validation.sh
#=========================================================================================================



source $CONDA_PREFIX/etc/profile.d/mamba.sh
source $CONDA_PREFIX/etc/profile.d/conda.sh

## Bnapus ZS11 setup

GFF=../validation/ZS11/zs11.genome.noScaf.gff
FASTA=../validation/ZS11/zs11.genome.fa
mamba activate prep
./prep --gff=$GFF --fasta=$FASTA --outdir=../validation/ZS11/data
mamba deactivate
## Bnapus ZS11 predictions

CUDA_VISIBLE_DEVICES=0 python nemo/validation.py ../validation/ZS11/data


## Athaliana setup

VAL=../validation/Athaliana
EQTL=${VAL}/data/eQTL
GATC=${VAL}/data/GATC
mkdir -p $EQTL
mkdir -p $GATC

GFF=../data/Athaliana/parent_data/TAIR10_GFF3_genes.gff

mamba activate prep
LIST=${VAL}/genes.eQTL.txt
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$EQTL"/"$STEM
  mkdir $OUTDIR
  echo $FASTA
  echo $OUTDIR
done | parallel -N 2 -j 10 ./prep --gff=$GFF --fasta={1} --outdir={2} --genelist=$LIST


# check if everythin was set up correctly
NUM_EXPECTED=$(wc -l $LIST | cut -f1 -d" ")
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$EQTL"/"$STEM
  if [[ $(wc -l $EQTL"/"$STEM/gene.id.key | cut -f1 -d" ") -ne $NUM_EXPECTED ]]; then
    rm -r $OUTDIR
    ./prep --gff=$GFF --fasta=$FASTA --outdir=$OUTDIR --genelist=$LIST
  fi
done

LIST=${VAL}/genes.GATC.txt
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$GATC"/"$STEM
  mkdir $OUTDIR
  echo $FASTA
  echo $OUTDIR
done | parallel -N 2 -j 10 ./prep --gff=$GFF --fasta={1} --outdir={2} --genelist=$LIST

# check if everythin was set up correctly
NUM_EXPECTED=$(wc -l $LIST | cut -f1 -d" ")
for FASTA in ${VAL}/*.fa
do
  STEM=$(basename $FASTA ".chr.fa")
  OUTDIR=$GATC"/"$STEM
  if [[ $(wc -l $GATC"/"$STEM/gene.id.key | cut -f1 -d" ") -ne $NUM_EXPECTED ]]; then
    rm -r $OUTDIR
    ./prep --gff=$GFF --fasta=$FASTA --outdir=$OUTDIR --genelist=$LIST
  fi
done
mamba deactivate

#### Athaliana predictions

export count=0; for file in ../validation/Athaliana/*.fa; do
  if [[ count -le 316 ]]; then
    echo $file >> ../validation/Athaliana/genomes_set1.txt
    let count++
  else
    echo $file >> ../validation/Athaliana/genomes_set2.txt
    let count++
  fi
done

modelfile=../model_weights/nemo/Athaliana/masked_graphpart/nemo_full.h5

mamba activate nemo
for fasta in $(cat ../validation/Athaliana/genomes_set1.txt)
do
    name=$(basename $fasta .chr.fa)
    dir=../validation/Athaliana/data/eQTL/$name
    CUDA_VISIBLE_DEVICES=0 python nemo/validation.py $dir $modelfile
    dir=../validation/Athaliana/data/GATC/$name
    CUDA_VISIBLE_DEVICES=0 python nemo/validation.py $dir $modelfile
done & pid1=$!

for fasta in $(cat ../validation/Athaliana/genomes_set2.txt)
do
    name=$(basename $fasta .chr.fa)
    dir=../validation/Athaliana/data/eQTL/$name
    CUDA_VISIBLE_DEVICES=1 python nemo/validation.py $dir $modelfile
    dir=../validation/Athaliana/data/GATC/$name
    CUDA_VISIBLE_DEVICES=1 python nemo/validation.py $dir $modelfile
done & pid2=$!

wait $pid1 $pid2
mamba deactivate


## Athaliana aggregate
OUTFILE=../validation/Athaliana/validation.eQTL.tsv
echo -e "accession_id\tgene_id\tpredicted_expression" > $OUTFILE
for DIR in ../validation/Athaliana/data/eQTL/*
do
  genotype=$(echo $DIR | rev | cut -f1 -d"/" | rev)
  FILE=$DIR"/predictions.txt"
  awk -v genotype="$genotype" 'BEGIN{OFS="\t"};NR>1{print genotype, $1, $2}' $FILE >> $OUTFILE
done

OUTFILE=../validation/Athaliana/validation.GATC.tsv
echo -e "accession_id\tgene_id\tpredicted_expression" > $OUTFILE
for DIR in ../validation/Athaliana/data/GATC/*
do
  genotype=$(echo $DIR | rev | cut -f1 -d"/" | rev)
  FILE=$DIR"/predictions.txt"
  awk -v genotype="$genotype" 'BEGIN{OFS="\t"};NR>1{print genotype, $1, $2}' $FILE >> $OUTFILE
done
