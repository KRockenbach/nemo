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
#title: setup_files.sh
#description: Sets-up and scales data partitions for training
#author: Kevin Rockenbach, Agnieszka Golicz
#email: kevin.rockenbach@ag.uni-giessen.de, agnieszka.golicz@agrar.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage:
#      - for nemo:           bash setup_files.sh <MASKED_GENOME_FASTA> <GENOME_FAI> <GENOME_GFF> <EXPRESSION_DATA> <DATADIR>
#      - for Xpresso + nemo: bash setup_files.sh <GENOME_FASTA> <MASKED_GENOME_FASTA> <ONLY_CDS_GENOME_FASTA> <GENOME_FAI> <GENOME_GFF> <EXPRESSION_DATA> <DATADIR> <CDS_FASTA>
#notes: intended to be run from within run_setup.sh
#       when running directly from command line, first activate prep environment
#=========================================================================================================

N_ARGS=$(echo $1 $2 $3 $4 $5 $6 $7 $8 | awk 'BEGIN{FS=" "}; {print NF}')

if [[ $N_ARGS -eq 8 ]]
then
  GENOME_FASTA=$1
  MASKED_GENOME_FASTA=$2
  ONLY_CDS_GENOME_FASTA=$3
  GENOME_FAI=$4
  GENOME_GFF=$5
  EXPRESSION=$6
  DATADIR=$7
  DERIVED=$DATADIR"/derived_data"
  CDS_FASTA=$8
  ORGANISM=$(echo $DATADIR | cut -f 3 -d "/")
else
  MASKED_GENOME_FASTA=$1
  GENOME_FAI=$2
  GENOME_GFF=$3
  EXPRESSION=$4
  DATADIR=$5
  DERIVED=$DATADIR"/derived_data"
  ORGANISM=$(echo $DATADIR | cut -f 3 -d "/")
fi


if [[ $N_ARGS -eq 8 ]]
then
  ## create a shuffled version of graphpart results ($DERIVED/graphpart_result_shuff.csv)
  python nemo/preprocessing/shuffle_graphpart_clusters.py $DERIVED/graphpart_result.csv
fi

###get start and end of interval around TSS and and TTS
###and store in promoters.bed, teminators.bed respectively

###genes where interval extends beyond end of chromosome are kept,
###to later extend intervals with Ns

###annotations for UTRs are redirected to noScaf.utr.gff

echo "Defining sequence interval bed files..."
##################################################
PROM_BED=$DERIVED"/promoters.bed"
TERM_BED=$DERIVED"/terminators.bed"
echo "Saving interval coordinates under $PROM_BED and $TERM_BED"
# $1 = outside interval; $2 = inside interval
python3 nemo/preprocessing/get_promoter_terminator.py 15000 5000 \
    $GENOME_FAI $GENOME_GFF mRNA $DERIVED
##################################################

if [[ $N_ARGS -eq 8 ]]
then
  echo "Retreiving UTR annotations"
  ###get 3p (TP) and 5p (FP) UTR data and save in seperate files
  FP_GFF=$(echo $GENOME_GFF | sed 's/.gff/.utr.5p.gff/' | sed 's/parent_/derived_/')
  TP_GFF=$(echo $GENOME_GFF | sed 's/.gff/.utr.3p.gff/' | sed 's/parent_/derived_/')
  echo "Saving UTR annotations under $FP_GFF and $TP_GFF"
  grep -P "\tfive_prime_UTR\t" $GENOME_GFF > $FP_GFF
  grep -P "\tthree_prime_UTR\t" $GENOME_GFF > $TP_GFF

  echo "Retreiving gene name, gene length and number of exons"
  GENE_STATS=$(echo $GENOME_GFF | sed 's/.gff/.mrna.stats/' | sed 's/parent_/derived_/')
  echo "Saving under $GENE_STATS"
  awk -f nemo/preprocessing/get_gene_stats.awk $GENOME_GFF > $GENE_STATS

  echo "Converting UTR annotations to BED format"
  TP_BED=$(echo $TP_GFF | sed 's/.gff/.bed/')
  FP_BED=$(echo $FP_GFF | sed 's/.gff/.bed/')
  echo "Saving BED files under $TP_BED and $FP_BED"
  python3 nemo/preprocessing/gff2bed.py $TP_GFF > $TP_BED
  python3 nemo/preprocessing/gff2bed.py $FP_GFF > $FP_BED

  echo "Retreiving UTR sequences"
  # get fasta files of UTRs (reverse complement if needed)
  FP_FASTA=$(echo $FP_BED | sed 's/.bed/.fasta/')
  TP_FASTA=$(echo $TP_BED | sed 's/.bed/.fasta/')
  echo "Saving UTR sequences under $TP_FASTA and $FP_FASTA"
  bedtools getfasta -fi $GENOME_FASTA -bed $FP_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $FP_FASTA
  bedtools getfasta -fi $GENOME_FASTA -bed $TP_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $TP_FASTA


  echo "Retreiving sequence length and GC contents of UTRs"
  FP_GC=$(echo $FP_FASTA | sed 's/.fasta/.gc/')
  TP_GC=$(echo $TP_FASTA | sed 's/.fasta/.gc/')
  echo "Saving data under $FP_GC and $TP_GC"
  # promoter.bed only needed for gene names
  python3 nemo/preprocessing/get_len_gc.py $FP_FASTA $PROM_BED > $FP_GC
  python3 nemo/preprocessing/get_len_gc.py $TP_FASTA $PROM_BED > $TP_GC

  echo "Retreiving sequence length and GC contents of CDS"
  CDS_GC=$(echo $CDS_FASTA | sed 's/.fasta/.gc/')
  echo "Saving data under $CDS_GC"
  # promoters.bed only needed for gene names
  python3 nemo/preprocessing/get_len_gc.py $CDS_FASTA $PROM_BED > $CDS_GC
fi

echo "Retreiving promoter and terminator sequences"
# get fasta file of promoters (-15kb TSS +5kb) and terminators (-5kb TTS +15kb)
# get reverse complement if needed
# Corresponding bed files created by get_promoter_terminator.py
if [[ $N_ARGS -eq 8 ]]
then
  PROM_FASTA=$(echo $PROM_BED | sed 's/.bed/.fasta/')
  TERM_FASTA=$(echo $TERM_BED | sed 's/.bed/.fasta/')
  bedtools getfasta -fi $GENOME_FASTA -bed $PROM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $PROM_FASTA
  bedtools getfasta -fi $GENOME_FASTA -bed $TERM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $TERM_FASTA
fi

MASKED_PROM_FASTA=$(echo $PROM_BED | sed 's/.bed/.cds_masked.fasta/')
MASKED_TERM_FASTA=$(echo $TERM_BED | sed 's/.bed/.cds_masked.fasta/')
bedtools getfasta -fi $MASKED_GENOME_FASTA -bed $PROM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $MASKED_PROM_FASTA
bedtools getfasta -fi $MASKED_GENOME_FASTA -bed $TERM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $MASKED_TERM_FASTA

if [[ $N_ARGS -eq 8 ]]
then
  ONLY_CDS_PROM_FASTA=$(echo $PROM_BED | sed 's/.bed/.only_cds.fasta/')
  ONLY_CDS_TERM_FASTA=$(echo $TERM_BED | sed 's/.bed/.only_cds.fasta/')
  bedtools getfasta -fi $ONLY_CDS_GENOME_FASTA -bed $PROM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $ONLY_CDS_PROM_FASTA
  bedtools getfasta -fi $ONLY_CDS_GENOME_FASTA -bed $TERM_BED -name -s | sed 's/(-)//' | sed 's/(+)//' > $ONLY_CDS_TERM_FASTA

  echo "Merging data together"
  for MASK_TYPE in "masked" "clear" "only_cds"
  do
    if [[ "$MASK_TYPE" == "clear" ]]; then
      PROM=$PROM_FASTA
      TERM=$TERM_FASTA
    elif [[ "$MASK_TYPE" == "masked" ]]; then
      PROM=$MASKED_PROM_FASTA
      TERM=$MASKED_TERM_FASTA
    else
      PROM=$ONLY_CDS_PROM_FASTA
      TERM=$ONLY_CDS_TERM_FASTA
    fi
    echo -n "$DERIVED $FP_GC $TP_GC $CDS_GC $GENE_STATS $EXPRESSION $PROM $TERM $MASK_TYPE "
  done | xargs -n 9 -P 3 python3 nemo/preprocessing/merge_data.py
else
  PROM=$MASKED_PROM_FASTA
  TERM=$MASKED_TERM_FASTA
  python3 nemo/preprocessing/merge_data.py $DERIVED $EXPRESSION $PROM $TERM masked
fi


cat ${DERIVED}"/gene.id.key.masked" > ${DERIVED}"/gene.id.key"
rm ${DERIVED}/gene.id.key.*


echo "Compressing data"
## compress merged data
gzip ${DERIVED}/merged.data*


if [[ $N_ARGS -eq 8 ]]
then
  # parallelized setup of datafolds
  for MASK_TYPE in "masked" "clear" "only_cds";
  do
    for CLUSTER_TYPE in "graphpart" "random";
    do
      FILE=$DERIVED"/merged.data.${MASK_TYPE}.gz"
      OUTDIR=$DATADIR"/${MASK_TYPE}_${CLUSTER_TYPE}_fold_data/"
      if [[ $CLUSTER_TYPE = "graphpart" ]]
      then
        GP_CSV=$DERIVED"/graphpart_result.csv"
      else
        GP_CSV=$DERIVED"/graphpart_result_shuff.csv"
      fi
      echo -n "$FILE $OUTDIR ${DERIVED}/gene.id.key $GP_CSV $MASK_TYPE $ORGANISM "
    done
  done | xargs -n 6 -P 6 python3 nemo/preprocessing/setup_folds.py
  python3 nemo/preprocessing/setup_full.py ${DERIVED}/merged.data.masked.gz ${DATADIR}/masked_graphpart_fold_data/
else
  FILE=$DERIVED"/merged.data.masked.gz"
  OUTDIR=$DATADIR"/masked_graphpart_fold_data/"
  GP_CSV=$DERIVED"/graphpart_result.csv"
  python3 nemo/preprocessing/setup_folds.py $FILE $OUTDIR ${DERIVED}/gene.id.key $GP_CSV masked $ORGANISM
  python3 nemo/preprocessing/setup_full.py ${DERIVED}/merged.data.masked.gz ${DATADIR}/masked_graphpart_fold_data/
fi


