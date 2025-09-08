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
#title: get_expression_class_lists.R
#description: classifies genes into three categories based on median expression (low, medium, high)
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript get_expression_class_lists.R <directory containing predictions> <output directory>
#=========================================================================================================

source("nemo/utils/eval_utils.R")

args <- commandArgs(trailingOnly=TRUE)

predfolder <- args[1]

df_path <- paste(predfolder, "/concat_preds.tsv", sep="")
df <- read.table(df_path, sep="\t", header=T)

expr_quantiles <- quantile(df$Actual, probs=c(0.25, 0.75))


transcript_names <- df$ID
names_lst <- strsplit(transcript_names, split = ".", fixed=TRUE)
gene_names <- c()
for (i in 1:length(transcript_names)){
  gene_names[i] <- names_lst[[i]][1]
}


high_IDs <- gene_names[df$Actual>=expr_quantiles[2]]
medium_IDs <- gene_names[df$Actual>expr_quantiles[1] & df$Actual<expr_quantiles[2]]
low_IDs <- gene_names[df$Actual<=expr_quantiles[1]]


IDfolder <- args[2]
if (!dir.exists(IDfolder)) dir.create(IDfolder, recursive=T)

fname <- paste(IDfolder, "/high_expr_ids.lst", sep="")
write.table(high_IDs, file=fname, quote=F, row.names=F, col.names=F, sep="\t")
fname <- paste(IDfolder, "/medium_expr_ids.lst", sep="")
write.table(medium_IDs, file=fname, quote=F, row.names=F, col.names=F, sep="\t")
fname <- paste(IDfolder, "/low_expr_ids.lst", sep="")
write.table(low_IDs, file=fname, quote=F, row.names=F, col.names=F, sep="\t")
