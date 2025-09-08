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


source("./plot_utils.R")



organisms <- c("Bnapus", "Athaliana")

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")


fig_lab_cex=1.5


plot_scatter <- function(organism){

  
  initial <- strsplit(organism, split="")[[1]][1]
  path <- paste0("../results/nemo/",organism,"/masked_graphpart/nemo",initial,"_preds/concat_preds.tsv")
  df <- read.table(path, header=T, sep="\t")
  
  x <- df$Actual
  y <- df$Predicted
  
  err <- (y - x)
  abs_err <- abs(err)
  cutoff <- median(abs_err)
  
  # Gaussian KDE
  df$Density <- get_density(df$Actual, df$Predicted, n = 100, h = c(1, 1))
  # scale density to range (0,1)
  df$Density <- df$Density/max(df$Density)
  # scale to integer range [1,1000]
  df$Density <- 1 + round(df$Density*999, 0)
  # sort df by density, so that densest points get drawn last
  df <- df[order(df$Density),]
  
  # 2A
  grey_scale <- c()
  for (i in 1:1000){
    grey_scale[i] <- grey(0.8/((i+199)/200))
  }
  par(mar=c(3,3,0,1), mgp=c(1.8,0.5,0), tck=-0.03)
  plot(df$Actual, df$Predicted, col=grey_scale[df$Density], pch=19,
       xlab = expression("Observed Expression [log"[10]*"(TPM + 0.1)]"),
       ylab = expression("Predicted Expression [log"[10]*"(TPM + 0.1)]  "),
       cex.lab=lab_cex, cex.axis=ax_cex,
       ylim=c(-1,4), xlim=c(-1,4))
  
  if (organism == "Athaliana"){
    legend("topleft", pch=NA, legend=expression(italic("A. thaliana")), text.col="grey40", bty="n")
  } else {
    legend("topleft", pch=NA, legend=expression(italic("B. napus")), text.col="grey40", bty="n")
  }
  
  Q <- quantile(x, probs=c(0.25,0.75))
  q1 <- Q[1]
  q3 <- Q[2]
  buffer <- median(abs(err))
  polygon(x=c(q1, q1, q3, q3),
          y=c(q1-buffer, q1+buffer, q3+buffer, q3-buffer),
          col=transparent("purple",0.65), border="purple")
  polygon(x=c(min(x)-0.05, min(x)-0.05, q1, q1),
          y=c(min(x)-0.05-buffer, min(x)-0.05+buffer, q1+buffer, q1-buffer),
          col=transparent("blue",0.65), border="blue")
  polygon(x=c(q3, q3, max(x)+0.05, max(x)+0.05),
          y=c(q3-buffer, q3+buffer, max(x)+0.05+buffer, max(x)+0.05-buffer),
          col=transparent("red",0.65), border="red")
  
  
  
  dx <- density(df$Actual, cut=F)
  par(mar=c(0,3,1,1), xpd=F)

  boxplot(df$Actual, type="n", axes=F, ylab="", xlab="", horizontal=T, ylim=c(-1,4))
  polygon(x=c(min(df$Actual)-0.05, min(df$Actual)-0.05, q1, q1), y=c(0,2,2,0), col=transparent("blue", 0.8), border=NA)
  polygon(x=c(q1, q1, q3, q3), y=c(0,2,2,0), col=transparent("purple", 0.8), border=NA)
  polygon(x=c(q3, q3, max(df$Actual)+0.05, max(df$Actual)+0.05), y=c(0,2,2,0), col=transparent("red", 0.8), border=NA)
  boxplot(df$Actual, axes=F, ylab="", xlab="", add=T, horizontal=T, col=NA)
  if (organism == "Bnapus"){
    fig_label(expression(bold("A")), cex=fig_lab_cex) 
  } else {
    fig_label(expression(bold("C")), cex=fig_lab_cex)
  }
}

plot_importance <-function(organism){
  initial <- strsplit(organism, split="")[[1]][1]
  p_path <- paste0("../results/nemo/",organism,"/masked_graphpart/attribs/promoter_expression_importance.tsv")
  t_path <- paste0("../results/nemo/",organism,"/masked_graphpart/attribs/terminator_expression_importance.tsv")
  pdf <- read.table(p_path, header=T, sep="\t")
  tdf <- read.table(t_path, header=T, sep="\t")
  for(c in c("high_expr_ids", "low_expr_ids", "medium_expr_ids")){
    pdf[,c] <- sgolayfilt(pdf[,c], p = 3, n = 111)
    tdf[,c] <- sgolayfilt(tdf[,c], p = 3, n = 111)
  }
  par(mar=c(4,4.2,2,0), xpd=F)
  plot(y=rep(0,6201),
       x=c(-5000:1200), type="n",
       xlab="", ylab="Importance", 
       axes=F, main="Promoter", cex.lab=lab_cex, cex.main=main_cex,
       ylim=c(0, max(c(pdf$high_expr_ids, pdf$low_expr_ids, pdf$medium_expr_ids,
                       tdf$high_expr_ids, tdf$low_expr_ids, tdf$medium_expr_ids))))
  if (organism == "Bnapus"){
    legend("topleft", pch=NA, legend=expression(italic("B. napus")), text.col="grey40", bty="n")
    fig_label(expression(bold("B")), cex=fig_lab_cex) 
  } else {
    legend("topleft", pch=NA, legend=expression(italic("A. thaliana")), text.col="grey40", bty="n")
    fig_label(expression(bold("D")), cex=fig_lab_cex)
  }
  axis(side=1, cex.axis=ax_cex, las=2)
  title(xlab="Position Relative to TSS", line=2.5)
  axis(side=2, cex.axis=ax_cex)
  lines(x=pdf$position, pdf$high_expr_ids, col="red", lwd=2)
  lines(x=pdf$position, pdf$medium_expr_ids, col="purple", lwd=2)
  lines(x=pdf$position, pdf$low_expr_ids, col="blue", lwd=2)
  lines(x=c(0,0), y=c(0,0.015), lwd=2, lty=2, col=transparent("grey30",0.5))
  
  par(mar=c(4,2,2,2.2))
  plot(y=rep(0,6201),
       x=c(-1200:5000), type="n",
       xlab="", ylab="", 
       cex.lab=lab_cex, axes=F, main="Terminator", cex.main=main_cex,
       ylim=c(0,max(c(pdf$high_expr_ids, pdf$low_expr_ids, pdf$medium_expr_ids,
                      tdf$high_expr_ids, tdf$low_expr_ids, tdf$medium_expr_ids))))
  axis(side=1, cex.axis=ax_cex, las=2)
  title(xlab="Position Relative to TTS", line=2.5)
  lines(x=tdf$position, tdf$high_expr_ids, col="red", lwd=2)
  lines(x=tdf$position, tdf$medium_expr_ids, col="purple", lwd=2)
  lines(x=tdf$position, tdf$low_expr_ids, col="blue", lwd=2)
  lines(x=c(0,0), y=c(0,0.015), lwd=2, lty=2, col=transparent("grey30",0.5))
} 



png("fig_3.png", height=14, width=17, res=1200, units="cm")
layout(matrix(c(2,3,4,
                1,3,4,
                6,7,8,
                5,7,8), ncol=3, byrow=T),
       heights=c(0.3,1,0.3,1),
       widths=c(1,1,1))

for (o in 1:2){
    organism <- organisms[o]
    seqs <- c("promoter", "terminator")
  

    lab_cex=1
    ax_cex=0.8
    main_cex=0.9
    
    plot_scatter(organism)
    plot_importance(organism)
    
}
dev.off()

