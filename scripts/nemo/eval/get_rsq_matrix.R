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


if(!require("stringr")){
    dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
    .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
    install.packages("stringr", repos = 'https://cloud.r-project.org')
}

library(stringr)

args <- commandArgs(trailingOnly=TRUE)

root <- ".."
model <- args[1]
organism <- args[2]
masking <- args[3]
partitioning <- args[4]
model_org <- args[5]


if (partitioning == "graphpart" | partitioning == "random"){
  all = TRUE
} else {
  all = FALSE
}
model_org_initial <- toupper(substr(model_org, 1, 1))


preddir <- paste(root, "/results/", model, "/", organism, "/", masking, "_", partitioning, "/", model, model_org_initial, "_preds", sep="")

if (all) {
  test_col <- c()
  valid_col <- c()
  N_col <- c()
  rsq_col <- c()
  org_col <- c()
  for (t in 0:9){
    predpath <- paste(preddir, "/predictions.t_", as.character(t), ".txt", sep="")
    actpath <- paste(preddir, "/actual.t_", as.character(t), ".txt", sep="")
    pred_df <- read.table(predpath, header=TRUE, sep="\t", row.names=1)
    act_df <- read.table(actpath, header=TRUE, sep="\t", row.names=1)
    pred_df <- as.data.frame(pred_df[order(rownames(pred_df)), ])
    act_df <- as.data.frame(act_df[order(rownames(act_df)), ])
    test_col = c(test_col, t)
    valid_col <- c(valid_col, NA)
    N_col <- c(N_col, NA)
    r <- cor(act_df[,1], pred_df[,1], method="pearson")
    rsq_col <- c(rsq_col, r**2)
    org_col <- c(org_col, organism)
    if ((masking == "masked") & (organism == "Bnapus") & (partitioning == "graphpart")){
      predpath_DS <- str_replace(string=predpath, pattern="_preds", replacement="_DS_preds")
      actpath_DS <- str_replace(string=actpath, pattern="_preds", replacement="_DS_preds")
      pred_df <- read.table(predpath_DS, header=TRUE, sep="\t", row.names=1)
      act_df <- read.table(actpath_DS, header=TRUE, sep="\t", row.names=1)
      pred_df <- as.data.frame(pred_df[order(rownames(pred_df)), ])
      act_df <- as.data.frame(act_df[order(rownames(act_df)), ])
      test_col <- c(test_col, t)
      valid_col <- c(valid_col, NA)
      N_col <- c(N_col, NA)
      r <- cor(act_df[,1], pred_df[,1], method="pearson")
      rsq_col <- c(rsq_col, r**2)
      org_col <- c(org_col, "BnapusDS")
    }
  }
  rsq_matrix <- cbind(test_col, valid_col, N_col, rsq_col)
  len <- length(rsq_col)
  rsq_matrix <- cbind(rep(model, len), org_col, rep(masking, len), rep(partitioning, len), rep(model_org, len), rsq_matrix)
} else {
  t <- 0
  test_col <- c()
  valid_col <- c()
  N_col <- c()
  rsq_col <- c()
  model_col <- c()
  if (model == "nemo"){
    validation <- c(as.character(1), "None")
  } else {
    validation <- c(as.character(1))
  }
  for (v in validation){
    for (N in 0:9){
      N <- as.character(N)
      # use stringr::str_replace(string, pattern, replacement) to adjust filenames conditionally
      # create N column and valid_fold column in df, save NA where not applicable
      predpath <- paste(preddir, "/predictions.t_", as.character(t), ".txt", sep="")
      actpath <- paste(preddir, "/actual.t_", as.character(t), ".txt", sep="")
      if (v != "None"){
        predpath <- str_replace(string=predpath, pattern=".txt", replacement=paste("_v_", v, ".txt", sep=""))
      }
      predpath <- str_replace(string=predpath, pattern=".txt", replacement=paste("_n_", N, ".txt", sep=""))
      pred_df <- read.table(predpath, header=TRUE, sep="\t", row.names=1)
      act_df <- read.table(actpath, header=TRUE, sep="\t", row.names=1)
      pred_df <- as.data.frame(pred_df[order(rownames(pred_df)), ])
      act_df <- as.data.frame(act_df[order(rownames(act_df)), ])
      test_col <- c(test_col, t)
      if (v == "None"){
        valid_col <- c(valid_col, NA)
      } else {
        valid_col <- c(valid_col, v)
      }
      N_col <- c(N_col, N) 
      r <- cor(act_df[,1], pred_df[,1], method="pearson")
      rsq_col <- c(rsq_col, r**2)
      model_col <- c(model_col, model)
      if (model == "nemo" & v != "None"){
        for (seq in c("prom", "term")){
          predpath <- paste(preddir, "/rand_", seq ,"_predictions.t_", as.character(t), "_v_", as.character(v), "_n_", as.character(N), ".txt", sep="")
          actpath <- paste(preddir, "/actual.t_", as.character(t), ".txt", sep="")
          pred_df <- read.table(predpath, header=TRUE, sep="\t", row.names=1)
          act_df <- read.table(actpath, header=TRUE, sep="\t", row.names=1)
          pred_df <- as.data.frame(pred_df[order(rownames(pred_df)), ])
          act_df <- as.data.frame(act_df[order(rownames(act_df)), ])
          test_col <- c(test_col, t)
          valid_col <- c(valid_col, v)
          N_col <- c(N_col, N)
          r <- cor(act_df[,1], pred_df[,1], method="pearson")
          rsq_col <- c(rsq_col, r**2)
          model_col <- c(model_col, paste(model, "_rand_", seq, sep=""))
        }
      }
    }
  }
  rsq_matrix <- cbind(as.integer(test_col), as.integer(valid_col), as.integer(N_col), as.numeric(rsq_col))
  len <- length(model_col)
  rsq_matrix <- cbind(model_col, rep(organism, len), rep(masking, len), rep(partitioning, len), rep(model_org, len), rsq_matrix)
}


colnames(rsq_matrix) <- c("model", "test_organism", "masking", "partitioning", "train_organism", "test_fold", "valid_fold", "N", "rsq")

out_file <- paste(root, "/results/rsq_df.tsv", sep="")
rsq_df <- as.matrix(read.table(out_file, header=TRUE, sep="\t"))
print(ncol(rsq_df))
print(ncol(rsq_matrix))
rsq_matrix <- rbind(rsq_df, rsq_matrix)
write.table(as.data.frame(rsq_matrix), file=out_file, row.names=F, col.names=T, quote=F, sep="\t")
