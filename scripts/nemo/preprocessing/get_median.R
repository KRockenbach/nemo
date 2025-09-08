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
#title: get_median.R
#description: calculated median expression across samples
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript get_median.R <expression_matrix> <output>
#=========================================================================================================

if(!require(dplyr)){
    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
    .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
    install.packages("dplyr", repos='http://cran.us.r-project.org')
}
library("dplyr")
     
args = commandArgs(trailingOnly=TRUE)

TPM_df <- read.table(args[1], header=T, comment.char="#", row.names=1, sep="\t", quote="")
# remove non-numeric columns
TPM_df <- select_if(TPM_df, is.numeric)
# replace NA with 0
for (c in 1:ncol(TPM_df)){
  TPM_df[is.na(TPM_df[,c]),c] <- 0
}


transcript <- rownames(TPM_df)
get_median <- function(x){
    return(median(x, na.rm=T))
}
median_expression <- apply(TPM_df, MARGIN=1, FUN=get_median)

OUT_df <- data.frame(transcript, median_expression)
write.table(OUT_df, file=args[2], 
            sep="\t", row.names=F, col.names=T, quote=F)

OUT_name <- sub(".tsv$", ".clean.tsv", args[1])
TPM_df <- as.data.frame(cbind(transcript,as.matrix(TPM_df)))
write.table(TPM_df, file=OUT_name,
            sep="\t", row.names=F, col.names=T, quote=F)
