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
#title: motif_insertion_pipeline.sh
#description: Inserts motifs of the repsective TF families into 100 random promoter or terminator sequences,
#             generates a data frame containing info on the gene, the motif and the location of insertion,
#             as well as the baseline prediction and new prediction after insertion
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: bash motif_insertion_pipeline.sh
#notes: commented links are temporary links and will not work (see filtering criteria)
#=========================================================================================================>

source ${CONDA_PREFIX}/etc/profile.d/mamba.sh
source ${CONDA_PREFIX}/etc/profile.d/conda.sh

mamba activate nemo

## JASPAR db was queried based on species = Arabidopsis thaliana / Arabidopsis lyrata subsp. lyrata, version = latest, collection = core, taxon = plants
## Some "families" are classified as "classes" under JASPAR (TCP, Homeobox, AP2EREBP, HSF, bHLH, bZIP, MADS)
## Note that the Jaspar database inconsistently categorized MYB(-related) factors as MYB(-related) or Myb(-related), this was accounted for
## For Arabidopsis lyatra there are only 5 profiles, all 5 are TCP

### Family: C4-GATA-related ---> 12 profiles
#wget http://jaspar.elixir.no/temp/20250811104538_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/C2C2gata_pfm.txt
### Family: GARP_G2-like ---> 20 profiles
#wget http://jaspar.elixir.no/temp/20250811104805_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/G2like_pfm.txt
### Family: NAC ---> 53 profiles
#wget http://jaspar.elixir.no/temp/20250811105007_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/NAC_pfm.txt
### Family: Trihelix ---> 12 profiles
#wget http://jaspar.elixir.no/temp/20250811105133_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/Trihelix_pfm.txt
### Class: basic helix-loop-helix factors (bHLH) ---> 40 profiles
#wget http://jaspar.elixir.no/temp/20250811105316_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/bHLH_pfm.txt
### Class: MADS box factors ---> 18 profiles
#wget http://jaspar.elixir.no/temp/20250811105601_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/MADS_pfm.txt
### Family: MYB ---> 57 profiles
#wget http://jaspar.elixir.no/temp/20250811105752_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/MYB_pfm.txt.subset1
### Family: Myb ---> 8 profiles
#wget http://jaspar.elixir.no/temp/20250811105920_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/MYB_pfm.txt.subset2
### concatenation ---> 65 profiles
#cat ../data/motifs/MYB_pfm.txt.subset1 ../data/motifs/MYB_pfm.txt.subset2 > ../data/motifs/MYB_pfm.txt
### Class: Basic leucine zipper factors (bZIP) ---> 38 profiles
#wget http://jaspar.elixir.no/temp/20250811110219_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/bZIP_pfm.txt
### Class: TCP
### Species: A. thaliana ---> 19 profiles
#wget http://jaspar.elixir.no/temp/20250811110528_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/TCP_pfm.txt.subset1
### Species: A. lyatra ---> 5 profiles
#wget http://jaspar.elixir.no/temp/20250811110711_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/TCP_pfm.txt.subset2
### concatenation ---> 24 profiles
#cat ../data/motifs/TCP_pfm.txt.subset1 ../data/motifs/TCP_pfm.txt.subset2 > ../data/motifs/TCP_pfm.txt
### Family: WRKY ---> 46 profiles
#wget http://jaspar.elixir.no/temp/20250811111025_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/WRKY_pfm.txt
### Family: DOF ---> 25 profiles
#wget http://jaspar.elixir.no/temp/20250811111447_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/C2C2dof_pfm.txt
### Class: Homeo domain factors ---> 31 profiles
#wget http://jaspar.elixir.no/temp/20250811111614_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/Homeobox_pfm.txt
### Family: MYB-related ---> 27 profiles
#wget http://jaspar.elixir.no/temp/20250811111951_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/MYBrelated_pfm.txt.subset1
### Family: Myb-related ---> 5 profiles
#wget http://jaspar.elixir.no/temp/20250811112108_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/MYBrelated_pfm.txt.subset2
### concatenation ---> 32 profiles
#cat ../data/motifs/MYBrelated_pfm.txt.subset1 ../data/motifs/MYBrelated_pfm.txt.subset2 > ../data/motifs/MYBrelated_pfm.txt
### Class: AP2/EREBP ---> 86 profiles
#wget http://jaspar.elixir.no/temp/20250811112247_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/AP2EREBP_pfm.txt
### Class: Heat shock factors ---> 10 profiles
#wget http://jaspar.elixir.no/temp/20250811112405_JASPAR2024_combined_matrices_305443_pfm.txt -O ../data/motifs/HSF_pfm.txt


#for FAMILY in C2C2gata G2like NAC Trihelix bHLH MADS MYB bZIP TCP WRKY C2C2dof Homeobox MYBrelated AP2EREBP HSF
#do
#  mkdir -p "../data/motifs/${FAMILY}"
#  echo "../data/motifs/${FAMILY}" | awk -f nemo/mutation/separate_motifs.awk - "../data/motifs/${FAMILY}_pfm.txt"
#done

for FAMILY in C2C2gata G2like NAC Trihelix bHLH MADS MYB bZIP TCP WRKY C2C2dof Homeobox MYBrelated AP2EREBP HSF
do
  At_ID_PATH=$(ls ../data/Athaliana/parent_data/TF_ids/* | grep "${FAMILY}\.")
  Bn_ID_PATH=$(ls ../data/Bnapus/parent_data/TF_ids/* | grep "${FAMILY}\.")
  CUDA_VISIBLE_DEVICES=0 python -m nemo.mutation.mutate_sequences -o Athaliana -a $At_ID_PATH -f $FAMILY --m "nemo" & pid1=$!
  CUDA_VISIBLE_DEVICES=1 python -m nemo.mutation.mutate_sequences -o Bnapus -a $Bn_ID_PATH -f $FAMILY -m "nemo" & pid2=$!
  wait $pid1 $pid2
done


mamba deactivate
