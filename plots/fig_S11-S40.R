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


source("plot_utils.R")

organisms <- c("Bnapus", "Athaliana")

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")



for (fam in families){
  full_fam <- families_full[which(families==fam)]
  for (o in 1:2){
    organism <- organisms[o]
    TFs <- as.vector(read.table(paste0("../results/nemo/", organism, 
                                       "/masked_graphpart/", fam, "_insertion/TF_names.lst"),
                                header=F, sep="\t")[,1])
    
    max_diffs <- c()
    for (t in 1:length(TFs)){
      TF <- TFs[t]
      seq <- "promoter"
      medium <- paste0("../results/nemo/", organism, "/masked_graphpart/", 
                       fam, "_insertion/", TF, "_medium_exp.tsv")
      medium_df <- read.table(medium, header=T, sep="\t")
      medium_df$prom_insert_idx <- (medium_df$prom_insert_idx - 5000)
      medium_df$term_insert_idx <- (medium_df$term_insert_idx - 1200)
      
      if (seq == "promoter"){
        medium_resid <- (medium_df$prom_mut_preds - medium_df$baseline)
      } else {
        medium_resid <- (medium_df$term_mut_preds - medium_df$baseline)
      }
      
      y <- c()
      x <- c()
      idx <- 1
      
      for (i in -1000:999){
        x[idx] <- i
        y[idx] <- median(medium_resid[medium_df$prom_insert_idx==i], na.rm=T)
        idx <- idx + 1
      }
      
      na_filt <- (!is.na(y))
      y <- y[na_filt]
      x <- x[na_filt]
      filtered <- sgolayfilt(y, p = 1, n = 75)
      
      max_diffs <- c(max_diffs, max(abs(filtered),na.rm=T))
    }
    names(max_diffs) <- TFs
    print(max_diffs)
    max_diffs <- sort(max_diffs, decreasing=T)
    if (length(max_diffs) >= 10){
      TFs <- names(max_diffs)[1:10]
    } else {
      TFs <- names(max_diffs)
    }
    print(TFs)
    # plot lineplots for TF subset with highest change in expression
    nrow=ceiling(length(TFs)/2)
    if (!dir.exists("fig_S11-S40")){
      dir.create("fig_S11-S40")
    }
    png(paste0("fig_S11-S40/",fam, "_", organism,"_insertion_lineplot.png"), width=17, height=(nrow*4), res=1200, units="cm")
    rows <- c(0,1,1)
    for (i in 1:nrow){
      offset <- (((i-1)*2)+3)
      rows <- c(rows,c(3,(offset+1),(offset+2)))
    }
    rows <- c(rows,c(0,2,2))
    layout(matrix(rows, byrow=T, ncol=3),
           widths=c(0.3,3,3), heights=c(0.3,rep(3,nrow), 0.3))
    
    par(mar=c(0,3,0,1))
    plot(x=c(1:100), y=c(1:100), type="n", axes=F, ylab="", xlab="")
    if (organism == "Athaliana"){
      text(x=50, y=50, labels=bquote(.(full_fam)~"("*italic("A. thaliana")*")"), col="grey30") 
    } else {
      text(x=50, y=50, labels=bquote(.(full_fam)~"("*italic("B. napus")*")"), col="grey30")
    }
    
    par(mar=c(0,3,0,1))
    plot(x=c(1:100), y=c(1:100), type="n", axes=F, ylab="", xlab="")
    text(x=50, y=50, labels=expression(bold("Position Relative to TSS")))
    
    
    par(mar=c(3,0,1,0))
    plot(x=c(1:100), y=c(1:100), type="n", axes=F, ylab="", xlab="")
    text(x=50, y=50, labels=expression(bold("Change in Expression Relative to Baseline")), srt=90)
    
    par(mar=c(3,3,1,1))
    
    for (t in 1:length(TFs)){
      TF <- TFs[t]
      seq <- "promoter"
      medium <- paste0("../results/nemo/", organism, "/masked_graphpart/", 
                       fam, "_insertion/", TF, "_medium_exp.tsv")
      medium_df <- read.table(medium, header=T, sep="\t")
      medium_df$prom_insert_idx <- (medium_df$prom_insert_idx - 5000)
      medium_df$term_insert_idx <- (medium_df$term_insert_idx - 1200)
      
      if (seq == "promoter"){
        medium_resid <- (medium_df$prom_mut_preds - medium_df$baseline)
      } else {
        medium_resid <- (medium_df$term_mut_preds - medium_df$baseline)
      }
      
      y <- c()
      yq1 <-c()
      yq3 <-c()
      x <- c()
      idx <- 1
      
      for (i in -1000:999){
        x[idx] <- i
        y[idx] <- median(medium_resid[medium_df$prom_insert_idx==i], na.rm=T)
        yq3[idx] <- quantile(medium_resid[medium_df$prom_insert_idx==i], probs=0.75, na.rm=T)
        yq1[idx] <- quantile(medium_resid[medium_df$prom_insert_idx==i], probs=0.25, na.rm=T)
        idx <- idx + 1
      }
      # filter NAs out
      na_filt <- (!is.na(y))
      y <- y[na_filt]
      yq1 <- yq1[na_filt]
      yq3 <- yq3[na_filt]
      x <- x[na_filt]
      na_filt <- (!is.na(yq1))
      y <- y[na_filt]
      yq1 <- yq1[na_filt]
      yq3 <- yq3[na_filt]
      x <- x[na_filt]
      na_filt <- (!is.na(yq3))
      y <- y[na_filt]
      yq1 <- yq1[na_filt]
      yq3 <- yq3[na_filt]
      x <- x[na_filt]
      filtered <- sgolayfilt(y, p = 1, n = 75)
      fq1 <- sgolayfilt(yq1, p = 1, n = 75)
      fq3 <- sgolayfilt(yq3, p = 1, n = 75)
      plot(x, filtered, type="l", ylim=c(min(fq1), max(fq3)), axes=F,
            ylab="",
            xlab="")
      polygon(x=c(x,rev(x)), y=c(fq3,rev(fq1)), col=transparent("blue", 0.9), border=NA)
      lines(x, filtered, lwd=2)
      abline(h=0, lty=2, col=transparent("grey40",0.5))
      if (max(filtered) <= 0.02){
        legend("bottomleft", pch=NA, legend=TF, bty="n", text.col="grey20") 
      } else {
        legend("topleft", pch=NA, legend=TF, bty="n", text.col="grey20")  
      }

      axis(2)
      axis(1)
    }
    
    dev.off()
    
  }  
}  


