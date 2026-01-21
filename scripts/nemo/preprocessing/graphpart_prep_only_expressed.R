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
#title: graphpart_prep_only_expressed.R
#description: Prepares data frame to be converted into graphpart input fasta file. Also outputs data frame of expression quartiles for all expressed genes
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2026-01-08
#version: 1.0.1
#usage: Rscript graphpart_prep_only_expressed.R <TSV of names and expression quartiles> <TSV of names and CDS sequence> <TSV of homoeolog pairs (optional)> <graphpart DF output> <expression quartile output> <species>
#=========================================================================================================


args = commandArgs(trailingOnly=TRUE)
#1 expression quartiles (.tsv)
#2 input "data/CDS_seqs.tsv"
#3 input "data/high_id_homoeolog_pairs.tsv" (optional)
#4 output "data/graphpart_df.tsv"
#5 output "data/only_expressed_quartiles_TPM.tsv"
#6 input organism

QUARTILE_TPM <- read.table(args[1], sep="\t", header=T)
MAX_TPM <- QUARTILE_TPM$max_expression
# exclude non-expressed genes (likely pseudogenes)
# they will not be used for model training and evaluation

only_expressed <- QUARTILE_TPM[MAX_TPM > 0,]
if (length(args)==6){
    # exclude high-ID homoeologs
    high_ID <- read.table(args[3], sep="\t", header=F)
    # concat together homoeologs from A and C subgenomes
    high_ID <- c(as.character(high_ID[,1]), as.character(high_ID[,2]))
    # exclude high ID homoeologs
    only_expressed <- only_expressed[!(only_expressed[,1] %in% high_ID),]
    gp_df_out <- args[4]
    tpm_out <- args[5]
    organism <- args[6] 
} else {
    gp_df_out <- args[3] 
    tpm_out <- args[4]
    organism <- args[5]
}

# 25% and 75% quantiles of expressed genes are boundaries for low, mid high expression classes
quant <- quantile(x=only_expressed$median_expression, probs=c(0.25,0.75)) 
class <- rep(paste(organism, "Medium", sep=""),nrow(only_expressed))
class_low <- only_expressed$median_expression <= quant[1]
class_mid <- (only_expressed$median_expression > quant[1] & only_expressed$median_expression < quant[2])
class_high <- only_expressed$median_expression >= quant[2]

class[class_low] <- paste(organism, "Low", sep="")
class[class_high] <- paste(organism, "High", sep="")
names(class) <- only_expressed$transcript

cds_seqs <- read.table(args[2], sep="\t", header=F)
colnames(cds_seqs) <- c('names', 'seqs')
names <- names(class)
class_df <- data.frame("names"=names, "class"=class)
merged_df <- merge(cds_seqs, class_df, by="names")

# since non-expressed and high-ID genes were already filtered out
# all remaining genes get priority 1
priority <- rep(1, nrow(merged_df))

cds_df <- cbind(merged_df, priority)
colnames(cds_df) <- c("names", "seqs", "class", "priority")
write.table(cds_df, 
            file=gp_df_out, 
            quote = F, sep="\t", dec=".", row.names=F, col.names = T)
# output median expression for genes categorized as expressed
# max expression > 0, median might still be 0
write.table(only_expressed,
            file=tpm_out,
            quote = F, sep="\t", dec=".", row.names=F, col.names = T)
