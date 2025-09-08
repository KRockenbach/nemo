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
#title: prep_Bnapus.sh
#description: Initial reformatting, standardization and preprocessing and partitioning of B. napus data
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: bash prep_Bnapus.sh
#notes: intended to be called from within run_setup.sh
#       when running directly from command line, first activate prep environment
#=========================================================================================================

ROOT=".."
DATADIR=${ROOT}"/data/Bnapus"
DERIVED=$DATADIR"/derived_data"
LOGS=$DATADIR"/logs"
PARENT=$DATADIR"/parent_data"
GENOME_FASTA=$PARENT"/Express617_v1.softmask.fa"
GENOME_GFF=$PARENT"/Express617_v1_Helixer.gff3"
SAMPLE_INFO=$PARENT"/sample_info_Express.txt"
EXPRESSION=$PARENT"/EAGLE_TPM_matrix.tsv"

rm -rf $DERIVED
mkdir -p $DERIVED
mkdir -p $LOGS

cp $GENOME_FASTA ${DERIVED}"/genome.fasta"
cp $GENOME_GFF ${DERIVED}"/genes.gff"
cp $EXPRESSION ${DERIVED}"/raw_expr_matrix.tsv"

GENOME_FASTA=${DERIVED}"/genome.fasta"
GENOME_GFF=${DERIVED}"/genes.gff"
EXPRESSION=${DERIVED}"/raw_expr_matrix.tsv"

echo "Preprocessing expression data..."
echo "Calculating averages across replicates..."
REP_AVRG=$DERIVED"/rep_avrg_matrix.tsv"
Rscript nemo/preprocessing/get_rep_avrg.R $SAMPLE_INFO $EXPRESSION $REP_AVRG
echo "Calculating medians expression across samples..."
MEDIANS=$DERIVED"/median_expr.tsv"
Rscript nemo/preprocessing/get_median.R $REP_AVRG $MEDIANS
rm $REP_AVRG

echo "Excluding unplaced contigs from annotation..."
NOSCAF_GFF=$DERIVED"/genes.filtered.gff"
grep "^chr" $GENOME_GFF > $NOSCAF_GFF

#echo "Reformating GFF file (standardizing UTR feature information)..."
# AGAT cannot output into same file as input
cat $NOSCAF_GFF > tmp_gff
rm $NOSCAF_GFF

agat_convert_sp_gxf2gxf.pl -g tmp_gff -o $NOSCAF_GFF
LOGNAME=$(basename $NOSCAF_GFF .gff)'.agat.log'
mv "tmp_gff.agat.log" $LOGS"/"$LOGNAME
rm tmp_gff

grep -P "\tmRNA\t" $NOSCAF_GFF | cut -f 9 | cut -d ";" -f 1 | cut -d "=" -f 2 > $DERIVED'/used_transcripts.lst'

# CDS sequences needed for graphpart partitioning and identification of high ID homoeologs
echo "Retreiving CDS sequences"
CDS_FASTA=${DERIVED}"/CDS.fasta"
echo "Saving CDS sequence under $CDS_FASTA"
gffread -x $CDS_FASTA -g $GENOME_FASTA $NOSCAF_GFF

########################################################################
# graphpart data partitioning (https://github.com/graph-part/graph-part)
########################################################################
# convert CDS fasta from multiline to singleline, so that sequences can be extracted from a single line following the sequence name
SINGLELINE_CDS_FASTA=$(echo $CDS_FASTA | sed 's/.fasta/.singleline.fasta/')
bash nemo/preprocessing/multiline_to_singleline_fasta.sh $CDS_FASTA $SINGLELINE_CDS_FASTA

########################
# identification of high ID genes
while read line
do
    A_hom=$(echo $line | cut -f1 -d " ")
    C_hom=$(echo $line | cut -f2 -d " ")
    A_head=$(echo '>'$A_hom)
    C_head=$(echo '>'$C_hom)

    grep -A1 -wF $A_head $SINGLELINE_CDS_FASTA > CDS.tmp.fa
    grep -A1 -wF $C_head $SINGLELINE_CDS_FASTA >> CDS.tmp.fa

    vsearch --allpairs_global CDS.tmp.fa --acceptall --blast6out CDS.tmp.out

    cat CDS.tmp.out >> $DERIVED"/CDS_align.out"

    rm CDS.tmp.fa
    rm CDS.tmp.out

# reciprocal best hits are taken from EAGLE-RC pipeline
done < $PARENT/A.vs.C.reciprocal_best
HIGH_ID=$DERIVED"/high_id_homoeolog_pairs.tsv"
awk 'BEGIN {OFS="\t"}; $3==100 {print $1,$2}' $DERIVED"/CDS_align.out" > $HIGH_ID
#######################

# get sequences
grep -v ">" $SINGLELINE_CDS_FASTA > seqs
# get names
grep ">" $SINGLELINE_CDS_FASTA | sed 's/>//g' > names
# create tsv of names and sequences
paste names seqs > $DERIVED"/CDS_seqs.tsv"
rm names seqs
# make dataframe with gene names, CDS sequences, classes and priorities --> outputs graphpart_df.tsv
# classes defined based on expression data
# secondary output is dataframe of median expression, containing only expressed genes (max TPM > 0)
OE_MEDIANS=$DERIVED"/only_expressed_median.tsv"
GP_DF_TSV=$DERIVED"/graphpart_df.tsv"
CDS_SEQS=$DERIVED"/CDS_seqs.tsv"
Rscript nemo/preprocessing/graphpart_prep_only_expressed.R $(echo $REP_AVRG | sed 's/.tsv/.clean.tsv/') $MEDIANS $CDS_SEQS $HIGH_ID $GP_DF_TSV $OE_MEDIANS "Bnapus"
# substitute all whitespace with comma --> turn into csv
GP_DF_CSV=$(echo $GP_DF_TSV | sed 's/.tsv/.csv/')
cat $GP_DF_TSV | tr -s '[:blank:]' ',' > $GP_DF_CSV
# create fasta file with custom headers, which graphpart takes as input
# csv_to_fasta.py is provided in the graphpart repository
python nemo/preprocessing/csv_to_fasta.py --file $GP_DF_CSV --identifier_col 0 --sequence_col 1 --label_col 2 --priority_col 3
# run graphpart for 10 partitions with threshold 0.3, recompute threshold 0.2, on 27 threads with triangular distance matrix for nucleotides
# 28 threads are available in total
echo "Grouping data into 10 evolutionarily independent partitions based on CDS similarity"
GP_FASTA=$(echo $GP_DF_CSV | sed 's/.csv/.fasta/')
graphpart mmseqs2needle -ff $GP_FASTA -th 0.3 -re 0.2 -pn priority -pa 10 -nt 27 -nu -tr
mv graphpart_result.csv graphpart_result_report.json $DERIVED"/graphpart_result.Bn.csv"
mv needleall.error $LOGS"/"

CDS_GFF=$DERIVED"/genes.filtered.CDS.gff"
grep -P "\tCDS\t" $NOSCAF_GFF > $CDS_GFF

echo "Masking Annotated CDS sequences from assembly"
# bedtools also accepts GFF files
MASKED_CDS_FASTA=$DERIVED"/genome.CDS_masked.fasta"
bedtools maskfasta -fi $GENOME_FASTA -bed $CDS_GFF -fo $MASKED_CDS_FASTA

echo "Masking non-CDS sequences from assembly"
GENOME_FAI=$GENOME_FASTA".fai"
cut -f 1,2 $GENOME_FAI | sort -k1,1 -k2,2n > $DERIVED/contig_lengths.txt
NON_CDS_BED=$DERIVED/non_CDS.bed
CDS_GFF_SORTED=$(echo $CDS_GFF | sed 's/.gff$/.sorted.gff/')
cat $CDS_GFF | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > $CDS_GFF_SORTED
bedtools complement -i $CDS_GFF_SORTED -g $DERIVED/contig_lengths.txt > $NON_CDS_BED
ONLY_CDS_FASTA=$DERIVED"/genome.only_CDS.fasta"
bedtools maskfasta -fi $GENOME_FASTA -bed $NON_CDS_BED -fo $ONLY_CDS_FASTA

Rscript nemo/preprocessing/clean_up_expression_matrix.R $EXPRESSION $DERIVED/expr_matrix.tsv

