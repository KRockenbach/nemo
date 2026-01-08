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
#title: prep_Athaliana.sh
#description: Initial reformatting, standardization and preprocessing of A. thaliana data
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2026-01-08
#version: 1.0.1
#usage: bash prep_Athaliana.sh
#notes: intended to be called from within run_setup.sh
#       when running directly from command line, first activate prep environment
#=========================================================================================================


ROOT=".."
DATADIR=$ROOT"/data/Athaliana"
DERIVED=$DATADIR"/derived_data"
LOGS=$DATADIR"/logs"
PARENT=$DATADIR"/parent_data"
GENOME_FASTA=$PARENT"/TAIR10_chr_all.fasta"
GENOME_GFF=$PARENT"/TAIR10_GFF3_genes.gff"
EXPRESSION=$PARENT"/E-CURD-1-query-results.tpms.tsv"

rm -rf $DERIVED
mkdir -p $DERIVED
mkdir -p $LOGS

cp $GENOME_FASTA ${DERIVED}"/genome.fasta"
cp $GENOME_GFF ${DERIVED}"/genes.gff"
cp $EXPRESSION ${DERIVED}"/raw_expr_matrix.tsv"

GENOME_FASTA=${DERIVED}"/genome.fasta"
GENOME_GFF=${DERIVED}"/genes.gff"
EXPRESSION=${DERIVED}"/raw_expr_matrix.tsv"


# fix column names of expression file
echo "Preparing expression data"
FIXED_HEADER_EXPRESSION=$(echo $EXPRESSION | sed 's/.tsv$/.fixed.tsv/')
head -n5 $EXPRESSION | tail -n1 | sed 's/ /_/g' | sed 's/,//g' > $FIXED_HEADER_EXPRESSION
tail -n+6 $EXPRESSION >> $FIXED_HEADER_EXPRESSION
# remove samples with low average spearman correlation
Rscript nemo/preprocessing/remove_low_corr_samples.R $FIXED_HEADER_EXPRESSION

echo "Calculating expression quartiles across samples..."
QUARTILES=$DERIVED"/quartile_expr.tsv"
# CONDENSED_EXPRESSION is matrix which only contains numerical columns except for gene name column (created by get_quartiles.R)
CONDENSED_EXPRESSION=$(echo $FIXED_HEADER_EXPRESSION | sed 's/.tsv$/.clean.tsv/')

USER=ubuntu
sudo usermod -a -G staff $USER
Rscript nemo/preprocessing/get_quartiles.R $FIXED_HEADER_EXPRESSION $QUARTILES


echo "Excluding chloroplast and mitochondrial contigs and comments from GFF"
ONLY_CHROM_GFF=$DERIVED"/genes.noScaf.gff"
awk '$1 != "ChrC" && $1 != "ChrM" && $1 !~ /^#/ {print $0}' $GENOME_GFF > $ONLY_CHROM_GFF


echo "Removing transcripts that don't contain a proper CDS"
ONLY_CODING_GFF=$(echo $ONLY_CHROM_GFF | sed 's/.gff$/.only_coding.gff/')
gffread -o $ONLY_CODING_GFF -g $GENOME_FASTA -CJVEF --keep-genes $ONLY_CHROM_GFF

echo "Selecting longest isoforms"
FILTERED_GFF=$(echo $ONLY_CODING_GFF | sed 's/.noScaf.only_coding.gff$/.filtered.gff/')
agat_sp_keep_longest_isoform.pl -g $ONLY_CODING_GFF -o $FILTERED_GFF
LOGNAME=genes.noScaf.only_coding.agat.log
mv $LOGNAME $LOGS"/"


# get list of used transcripts
grep -P "\tmRNA\t" $FILTERED_GFF | cut -f 9 | cut -d ";" -f 1 | cut -d "=" -f 2 > $DERIVED'/used_transcripts.lst'


SINGLELINE_GENOME_FASTA=$(echo $GENOME_FASTA | sed 's/.fasta$/.singleline.fasta/')
bash nemo/preprocessing/multiline_to_singleline_fasta.sh $GENOME_FASTA $SINGLELINE_GENOME_FASTA
cat $SINGLELINE_GENOME_FASTA > $GENOME_FASTA
rm $SINGLELINE_GENOME_FASTA

GENOME_FAI=$(echo $GENOME_FASTA | sed 's/\.fasta$/\.fasta.fai/')
samtools faidx -n 0 -o $GENOME_FAI $GENOME_FASTA

echo "Retreiving CDS sequences"
# Partitioning is done based on CDS sequences
CDS_LONG_ISO_FASTA=${DERIVED}"/CDS.fasta"
echo "Saving CDS sequence under $CDS_LONG_ISO_FASTA"
# used for graphpart and masking
gffread -x $CDS_LONG_ISO_FASTA -g $GENOME_FASTA $FILTERED_GFF
# Convert transcript names to gene names
cut -d "." -f 1 $CDS_LONG_ISO_FASTA > tmp
cat tmp > $CDS_LONG_ISO_FASTA
rm tmp
##########################


########################################################################
# graphpart data partitioning (https://github.com/graph-part/graph-part)
########################################################################
# convert CDS fasta from multiline to singleline, so that sequences can be extracted from a single line following the sequence name
SINGLELINE_CDS_LONG_ISO_FASTA=$(echo $CDS_LONG_ISO_FASTA | sed 's/.fasta/.singleline.fasta/')
bash nemo/preprocessing/multiline_to_singleline_fasta.sh $CDS_LONG_ISO_FASTA $SINGLELINE_CDS_LONG_ISO_FASTA

# get sequences
grep -v ">" $SINGLELINE_CDS_LONG_ISO_FASTA > seqs
# get names
grep ">" $SINGLELINE_CDS_LONG_ISO_FASTA | sed 's/>//g' | cut -d "." -f 1 > names
# create tsv of names and sequences
paste names seqs > $DERIVED"/CDS_seqs.tsv"
rm names seqs
# make dataframe with gene names, CDS sequences, classes and priorities --> outputs graphpart_df.tsv
# classes defined based on expression data
# secondary output is dataframe of expression quartiles, containing only expressed genes (max TPM > 0)
OE_QUARTILES=$DERIVED"/only_expressed_quartiles.tsv"
GP_DF_TSV=$DERIVED"/graphpart_df.tsv"
CDS_SEQS=$DERIVED"/CDS_seqs.tsv"
Rscript nemo/preprocessing/graphpart_prep_only_expressed.R $QUARTILES $CDS_SEQS $GP_DF_TSV $OE_QUARTILES "Athaliana"
# substitute all whitespace with comma --> turn into csv
GP_DF_CSV=$(echo $GP_DF_TSV | sed 's/.tsv/.csv/')
cat $GP_DF_TSV | tr -s '[:blank:]' ',' > $GP_DF_CSV
# create fasta file with custom headers, which graphpart takes as input
# csv_to_fasta.py is provided in the graphpart repository (https://github.com/graph-part/graph-part)
python nemo/preprocessing/csv_to_fasta.py --file $GP_DF_CSV --identifier_col 0 --sequence_col 1 --label_col 2 --priority_col 3

CDS_GFF=$(echo $FILTERED_GFF | sed 's/.gff$/.CDS.gff/')
MASKED_FASTA=$(echo $GENOME_FASTA | sed 's/\.fasta$/\.CDS_masked.fasta/')
awk '$3=="CDS" {print $0}' $FILTERED_GFF > $CDS_GFF

echo "Masking annotated CDS sequences from assembly"
# bedtools also accepts GFF files
bedtools maskfasta -fi $GENOME_FASTA -bed $CDS_GFF -fo $MASKED_FASTA


Rscript nemo/preprocessing/clean_up_expression_matrix.R $CONDENSED_EXPRESSION $DERIVED/expr_matrix.tsv


