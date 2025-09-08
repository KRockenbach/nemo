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
#title: clean_up_expression_matrix.R
#description: some basic clean up of the expression matrix
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript clean_up_expression_matrix.R <expression_matrix> <output>
#=========================================================================================================

if(!require(dplyr)){
    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
    .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
    install.packages("dplyr", repos='http://cran.us.r-project.org')
}
library("dplyr")

args = commandArgs(trailingOnly=TRUE)

expr_path <- args[1]
out_path <- args[2]

expr_matrix <- read.table(expr_path, header=T, sep="\t", comment.char="#", quote="", row.names=1)
# get rid of extra columns
expr_matrix <- select_if(expr_matrix, is.numeric)

get_max <- function(x){
    return(max(x, na.rm=T))
}
max_expr <- apply(expr_matrix, MARGIN = 1, FUN=get_max)

# exclude non-expressed genes
only_expressed <- expr_matrix[max_expr > 0,]

Gene <- rownames(only_expressed)
# make sure names are gene-level names
names_lst <- strsplit(Gene, split = ".", fixed=TRUE)
gene_names <- c()
for (i in 1:length(Gene)){
  gene_names[i] <- names_lst[[i]][1]
}
Gene <- gene_names

only_expressed <- as.data.frame(cbind(Gene, as.matrix(only_expressed))) 

write.table(only_expressed, file=out_path, quote=F, row.names=F, col.names=T, sep="\t")
