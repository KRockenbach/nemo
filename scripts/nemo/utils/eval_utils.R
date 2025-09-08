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
#title: eval_utils.R
#description: helper functions for categorizing genes based on expression, specificity and predictions
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: from within R script, import functions using: source("relative/path/to/eval_utils.R")
#=========================================================================================================



##++++## ESTIMATE PREDICTION ERRORS ##+++++##

get_unity_error <- function(df){
  x <- df$Actual
  y <- df$Predicted
  # signed error (overpredictions are larger than 0)
  errs <- (y - x)
  names(errs) <- df$ID
  return(errs)
}

get_quantile_error <- function(df){
  x <- df$Actual
  y <- df$Predicted
  qx <- quantile(x, probs=c(0.25,0.75))
  qy <- quantile(y, probs=c(0.25,0.75))
  
  # delta y / delta x
  m <- (qy[2]-qy[1])/(qx[2]-qx[1])
  print(paste("m = ", as.character(m), sep=""))
  # y=mx+b
  # b=y-mx
  b <- qy[1] - (m*qx[1])
  print(paste("b = ", as.character(b), sep=""))
  # signed error (overpredictions are larger than 0)
  errs <- (y - ((m*x)+b))
  names(errs) <- df$ID
  return(errs)
}

##++++## ID CONVERSION ##++++##

convert_IDs <- function(IDs, key_path){
  key_df <- read.table(key_path, # gene.id.key
                       sep="\t", header=F)
  key_df <- key_df[key_df[,2] %in% IDs,]
  gene_names <- key_df[match(IDs,key_df[,2]),1]
  return(gene_names)
}


convert_names <- function(names, key_path){
  key_df <- read.table(key_path, # gene.id.key
                       sep="\t", header=F)
  key_df <- key_df[key_df[,1] %in% names,]
  gene_keys <- key_df[match(names,key_df[,1]),2]
  return(gene_keys)
}

