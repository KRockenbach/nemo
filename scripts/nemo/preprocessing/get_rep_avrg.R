#!/usr/bin/Rscript --vanilla

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


#=========================================================================================================
#title: get_rep_avrg.R
#description: calculates average expression across replicates of the same sample
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript get_rep_avrg.R <sample_information_table> <expression_matrix> <output_filename>
#=========================================================================================================

args = commandArgs(trailingOnly=TRUE)

sample_info <- read.table(args[1], sep="\t", header=T)

TPM_matrix <- read.table(args[2], sep="\t", header=T)

# initialize with gene names
rep_avrg_matrix <- as.character(TPM_matrix[,1])

subtissues <- as.character(unique(sample_info$subtissue))
for (st in 1:length(subtissues)){
  subtiss <- subtissues[st]
  samples <- sample_info$sample_name[sample_info$subtissue == subtiss]
  sub_matrix <- as.matrix(TPM_matrix[,colnames(TPM_matrix) %in% samples])
  n_reps <- sum(colnames(TPM_matrix) %in% samples)
  print(paste("Averaging over ", 
              as.character(n_reps), 
              " replicates for subtissue ", 
              subtiss, 
              sep=""))
  subtiss_avrg <- apply(sub_matrix, MARGIN=1, FUN=mean)
  rep_avrg_matrix <- cbind(rep_avrg_matrix, subtiss_avrg)
}

colnames(rep_avrg_matrix) <- c("transcript", subtissues)

write.table(rep_avrg_matrix, file=args[3],
            col.names=T, row.names=F, quote=F, sep="\t")
