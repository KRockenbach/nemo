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
#title: run_setup.sh
#description: runs the entire preprocessing pipeline, partitions data, sets up training files and fits scalers
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2026-01-08
#version: 1.0.1
#usage: bash run_setup.sh
#=========================================================================================================

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate prep


ROOT=".."
mkdir -p ${ROOT}"/data/Athaliana" ${ROOT}"/data/Bnapus"

bash prep_Athaliana.sh & pid1=$!
bash prep_Bnapus.sh & pid2=$!

wait $pid1 $pid2


### run graphpart for 10 partitions with threshold 0.3, recompute threshold 0.2, on 27 threads with triangular distance matrix for nucleotides
##echo "Grouping data into 10 evolutionarily independent partitions based on CDS similarity"
##cat ../data/Bnapus/derived_data/graphpart_df.fasta ../data/Athaliana/derived_data/graphpart_df.fasta > tmp.fasta
##graphpart mmseqs2needle -ff tmp.fasta -th 0.3 -re 0.2 -pn priority -pa 10 -nt $(nproc) -nu -tr

At_GP_RESULT="${ROOT}/data/Athaliana/derived_data/graphpart_result.csv"
Bn_GP_RESULT="${ROOT}/data/Bnapus/derived_data/graphpart_result.csv"
##head -n 1 graphpart_result.csv > $At_GP_RESULT
##head -n 1 graphpart_result.csv > $Bn_GP_RESULT
##tail -n +2 graphpart_result.csv | grep -v "Bnapus" >> $At_GP_RESULT
##tail -n +2 graphpart_result.csv | grep "Bnapus" >> $Bn_GP_RESULT
##mv graphpart_result* needleall.error $ROOT"/data/"
##rm tmp.fasta
cp "${ROOT}/data/Athaliana/parent_data/graphpart_result.joint.csv" $At_GP_RESULT
cp "${ROOT}/data/Bnapus/parent_data/graphpart_result.joint.csv" $Bn_GP_RESULT
cp "${ROOT}/data/Bnapus/parent_data/graphpart_result.Bn.csv" "${ROOT}/data/Bnapus/derived/"

for ORGANISM in "Bnapus" "Athaliana"
do
  #7
  DATADIR="${ROOT}/data/${ORGANISM}"
  DERIVED=$DATADIR"/derived_data"
  #$1
  GENOME_FASTA=""
  #$3
  ONLY_CDS_GENOME_FASTA=""
  #$8
  CDS_FASTA=""
  if [[ $ORGANISM = "Bnapus" ]]
  then
    GENOME_FASTA=$DERIVED"/genome.fasta"
    ONLY_CDS_GENOME_FASTA=$DERIVED"/genome.only_CDS.fasta"
    CDS_FASTA=$DERIVED"/CDS.singleline.fasta"
  fi
  #$2
  MASKED_GENOME_FASTA=$DERIVED"/genome.CDS_masked.fasta"
  #$4
  GENOME_FAI=$DERIVED"/genome.fasta.fai"
  #$5
  GENOME_GFF=$DERIVED"/genes.filtered.gff"
  #$6
  EXPRESSION=$DERIVED"/only_expressed_quartiles.tsv"
  # space at the end is important
  bash setup_files.sh $GENOME_FASTA $MASKED_GENOME_FASTA $ONLY_CDS_GENOME_FASTA $GENOME_FAI $GENOME_GFF $EXPRESSION $DATADIR $CDS_FASTA
done

bash downsample_data.sh

mamba deactivate
