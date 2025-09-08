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
#title: over_under_pred.R
#description: classifies genes into three categories based on prediction (overpredicted, well-predicted, underpredicted)
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript over_under_pred.R <directory containing predictions> <output directoy> 
#=========================================================================================================

source("nemo/utils/eval_utils.R") #relative to bash script

args <- commandArgs(trailingOnly=TRUE)

predfolder <- args[1]
input_fname <- paste(predfolder, "/concat_preds.tsv", sep="")
df <- read.table(input_fname, sep="\t", header=T)


# prediction errors
err <- get_unity_error(df=df)
abs_err <- abs(err)
act_vals <- df$Actual
pred_vals <- df$Predicted
cutoff <- median(abs_err)


well_keys <- names(err)[abs_err<=cutoff]
over_keys <- names(err)[err>cutoff]
under_keys <- names(err)[err<(-cutoff)]


# create output folder
IDfolder <- args[2]
if (!dir.exists(IDfolder)) dir.create(IDfolder, recursive=T)

# write output
write.table(well_keys, file=paste(IDfolder, "/well-predicted.lst", sep=""), quote=F,
               row.names=F, col.names=F, sep="\t")
write.table(over_keys, file=paste(IDfolder, "/overpredicted.lst", sep=""), quote=F,
               row.names=F, col.names=F, sep="\t")
write.table(under_keys, file=paste(IDfolder, "/underpredicted.lst", sep=""), quote=F,
               row.names=F, col.names=F, sep="\t")
