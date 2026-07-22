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

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")


get_mutation_data <- function(organism, type="TFBS"){
  pdf <- data.frame(position=c(-500:499))
  tdf <- data.frame(position=c(-500:499))
  for (family in families){
    df <- read.table(paste0("../result_subset/nemo/",organism,"/masked_graphpart/motif_insertion/", family, "_medium.tsv"), header=T)
    df <- df[df$motif_type==type,]
    df$diff <- (df$mutated_pred - df$baseline)
    pdiff <- c()
    tdiff <- c()
    for (p in pdf$position){
      pdiff <- c(pdiff, median(df$diff[df$mutated_sequence_type=="promoter" & df$relative_idx==p]))
      tdiff <- c(tdiff, median(df$diff[df$mutated_sequence_type=="terminator" & df$relative_idx==p]))
    }
    pdf <- cbind(pdf, pdiff)
    colnames(pdf)[ncol(pdf)] <- family
    tdf <- cbind(tdf, tdiff)
    colnames(tdf)[ncol(tdf)] <- family
  }
  return(list(pdf,tdf))
}


png("fig_S28.png", height=20, width=17, units="cm", res=1200)
par(mfrow=c(2,1), mar=c(4.2,4.2,0.1,0.1))

for (organism in c("Bnapus", "Athaliana")){

  df <- get_mutation_data(organism, "TFBS")
  idf <- df[[1]]
      
  auc <- function(x){
    #x <- x[501:1501]
    area <- 0
    for (i in 1:length(x)){
      if (i == 1){
        prev <- x[i]
        next
      } else {
        current <- x[i]
        # is previous > current, third term becomes negative
        # negative values result in negarive area
        area <- area + prev + ((current-prev)/2)
      }
    }
    return(area)
  }
  
  idf <- idf[,-1]
  area <- apply(as.matrix(idf), MARGIN=2, FUN=auc)
  area <- area[order(names(area))]
  
  initial <- strsplit(organism, split="")[[1]][1]
  exp_path <- paste0("../result_subset/nemo/", organism, "/masked_graphpart/nemo", initial, "_preds/concat_preds.tsv")
  exp_df <- read.table(exp_path, header=T, sep="\t")
  
  families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
                "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
  upstream <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
                "TCP", "bZIP", "MYB", "MADS")
  both <- c("bHLH", "Trihelix", "NAC", "G2like")
  downstream <- c("HSF", "AP2EREBP", "C2C2gata")
  
  ppath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.promoter.tsv"
  tpath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.terminator.tsv"
  
  
  ###################################
  
  # load and transpose df
  load_and_transpose <- function(path){
    df <- read.table(path, header=F)
    col_vec <- df[,1]
    df <- df[,-1]
    df <- t(df)
    df <- as.data.frame(df)
    colnames(df) <- col_vec
    rownames(df) <- c(1:nrow(df))
    return(df)
  }
  pdf <- load_and_transpose(ppath)
  
  pdf <- pdf[,families]
  
  # scale data within families
  for (f in families){
    norm <- min_max.positive.scale(pdf[,f])
  }
  
  weighted_avg_pos <- c()
  for (f in families){
    weighted_avg_pos[f] <- sum(pdf[,f] * c(-500:500))/sum(pdf[,f])
  }
  
  
  ##############################
  exp <- c()
  index <- 1
  for (f in families){
    if (organism == "Bnapus"){
      ID_path <- paste0("../data_subset/Bnapus/parent_data/TF_ids/",f,".at.cov.ids")
    } else {
      ID_path <- paste0("../data_subset/Athaliana/parent_data/TF_ids/",f,".all.narrowPeak.n.intersect.ids")
    }
    IDs <- read.table(ID_path, sep="\t", header=F)[,1]
    for (i in 1:length(IDs)){
      IDs[i] <- strsplit(IDs[i], split=".", fixed=T)[[1]][1]
    }
    exp[index] <- mean(exp_df$Actual[exp_df$ID %in% IDs])
    index <- index + 1
  }
  
  names(exp) <- families
  
  
  
  cols <- c("brown","red","tomato", "tan1", "gold", "khaki", "lawngreen",
            "aquamarine", "seagreen", "cornflowerblue", "blue", "darkslategray", "plum", "purple", "magenta")

  names(cols) <- families
  
  plot(weighted_avg_pos, exp, col=cols, pch=19, cex=2, bty="n",
       ylab="Average Expression Level",
       xlab='')
  
  if (organism == "Athaliana"){
    title(xlab='"Center of Mass" of DAP-Seq Peak Coverage Relative to TSS')
  }
  
  
  legend("topleft", legend = families_full, pch=19, col=cols, ncol=2, cex=0.8)
  
  
  abline(lm(exp~weighted_avg_pos))
  ct <- cor.test(exp,weighted_avg_pos, method="pearson")
  r <- ct$estimate
  p <- ct$p.value
  
  if (organism == "Bnapus"){
    org_lab <- "B. napus"
  } else {
    org_lab <- "A. thaliana"
  }
  legend("bottomright", bty="n", legend=bquote(italic(.(org_lab))), text.col="grey40")
  Lines <- list(bquote(rho["p"]*" = "*.(round(r,3))),bquote("p-value = "*.(round(p,3))))
  legend("top", bty="n", legend=do.call(expression,Lines)) 

}


dev.off()

