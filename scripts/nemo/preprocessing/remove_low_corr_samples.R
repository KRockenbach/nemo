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
#title: remove_low_corr_samples.R
#description: correlated expression of samples based on 5000 most highly expressed genes, then removes samples with average spearman correlation below or equal to 0.1
#author: Kevin Rockenbach
#email: kevin.rockenbach@ag.uni-giessen.de
#date: 2025-08-28
#version: 1.0.0
#usage: Rscript remove_low_corr_samples.R <expression_matrix>
#=========================================================================================================


if(!require(dplyr)){
    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
    .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
    install.packages("dplyr", repos='http://cran.us.r-project.org')
}
library("dplyr")


args = commandArgs(trailingOnly=TRUE)
df <- read.table(args[1], header=T, sep="\t", comment.char="#", quote="")


Gene_ID <- df[,1] 
Gene_Name <- df[,2]
# remove non-numeric columns
df <- select_if(df, is.numeric)
# replace NA with 0
for (c in 1:ncol(df)){
  df[is.na(df[,c]),c] <- 0
}


rowsums <- apply(df, MARGIN=1, FUN=sum)
top5000 <- order(rowsums, decreasing=T)

sub <- df[top5000,]
sub <- sub[1:5000,]


avrg_cors <- c()


for(c in 1:ncol(sub)){
  cors <- c()
  for(r in 1:ncol(sub)){
    cors[r] <- cor(sub[,c], sub[,r], method="spearman")
  }
  avrg_cors[c] <- mean(cors)
} 

names(avrg_cors) <- colnames(df)
to_be_removed <- names(avrg_cors[avrg_cors <= 0.1])
print(to_be_removed)

print(length(colnames(df)))
df <- df[,which(!(colnames(df) %in% to_be_removed))]
print(length(colnames(df)))


df <- cbind(Gene_ID, Gene_Name, df)

write.table(df, file=args[1], row.names=F, col.names=T, quote=F, sep="\t")
