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




families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")




###### FIGURES S11-S16 ##### normalized importance

ppath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.promoter.tsv"
tpath <- "../data_subset/plot_data/DAPseq/DAPSeq.TF.terminator.tsv"

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
tdf <- load_and_transpose(tpath)
tdf <- tdf[c(nrow(tdf):1),] # positions of terminator df are reversed

pdf <- pdf[,families]
tdf <- tdf[,families]

for (f in families){
  R=1
  N=101
  P=1
  pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
  tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
}

dap_p <- pdf
dap_t <- tdf


get_importance <- function(organism){
  pdf <- read.table(paste0("../result_subset/nemo/",organism,
                           "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                    header=T)
  tdf <- read.table(paste0("../result_subset/nemo/",organism,
                           "/masked_graphpart/attribs/terminator_TF_importance.tsv"),
                    header=T)
  tdf$position <- tdf$position + 1 # shift position, because TTS is part of upstream sequence
  # trim to +- 500 bp
  pdf <- pdf[501:1501,]
  tdf <- tdf[501:1501,]
  pdf <- pdf[,c("position","mean", "sd", families)]
  tdf <- tdf[,c("position","mean", "sd", families)]
  
  # mean normalize
  for (f in families){
    concat <- c(pdf[,f], tdf[,f])
    concat <- mean_norm(concat) #!
    pdf[,f] <- concat[1:length(pdf[,f])]
    tdf[,f] <- concat[(length(pdf[,f])+1):length(concat)]
    
  }
  # mean normalize (across sequence) mean (across samples)
  concat_mean <- c(pdf$mean, tdf$mean)
  concat_mean <- mean_norm(concat_mean) #!
  pdf$mean <- concat_mean[1:length(pdf$mean)]
  tdf$mean <- concat_mean[(length(pdf$mean)+1):length(concat_mean)]
  
  
  # standardize
  for (f in families){
    pdf[,f] <- (pdf[,f] - pdf$mean)
    tdf[,f] <- (tdf[,f] - tdf$mean)
  }
  
  for (f in families){
    pdf[,f] <- (pdf[,f]/pdf$sd)
    tdf[,f] <- (tdf[,f]/tdf$sd)
  }
  
  # savitzky-golay filtering
  for (c in 2:ncol(pdf)){
    R=1
    N=101 #75
    P=1
    pdf[,c] <- sgolayfilt(pdf[,c], p = P, n = N)
    tdf[,c] <- sgolayfilt(tdf[,c], p = P, n = N)
  }
  
  return(list(pdf,tdf))
}

At_imp <- get_importance("Athaliana")
At_imp_p <- At_imp[[1]]
At_imp_t <- At_imp[[2]]
Bn_imp <- get_importance("Bnapus")
Bn_imp_p <- Bn_imp[[1]]
Bn_imp_t <- Bn_imp[[2]]




plot_lines <- function(imp_p, imp_t, dap_p, dap_t, family){
  # promoter
  par(mar=c(2,5,1,0.5))
  plot(c(-500,500),c(max(c(max(dap_p),max(dap_t))),min(c(min(dap_p),min(dap_t)))), 
       type="n", bty="n", axes=F, ylab="", xlab="")
  lines(x=c(-500:500), y=dap_p, col="red", lwd=2)
  axis(1)
  axis(2, col="red", col.axis="red")
  title(ylab="DAP-Seq Peak Coverage", col.lab="red")
  par(new=T)
  plot(c(-500,500),c(max(c(max(imp_p),max(imp_t))),min(c(min(imp_p),min(imp_t)))), 
       type="n", bty="n", axes=F, ylab="", xlab="")
  lines(x=c(-500:500), y=imp_p, col="black", lwd=2)
  legend("topleft", legend=family, cex=1, pch=NA)

  
  #terminator
  par(mar=c(2,0.5,1,5))
  plot(c(-500,500),c(max(c(max(dap_p),max(dap_t))),min(c(min(dap_p),min(dap_t)))), 
       type="n", bty="n", axes=F, ylab="", xlab="")
  lines(x=c(-500:500), y=dap_t, col="red", lwd=2)
  axis(1)
  par(new=T)
  plot(c(-500,500),c(max(c(max(imp_p),max(imp_t))),min(c(min(imp_p),min(imp_t)))), 
       type="n", bty="n", axes=F, ylab="", xlab="")
  lines(x=c(-500:500), y=imp_t, col="black", lwd=2)
  axis(4)
  mtext("Normalized Importance", side=4, line=3, cex=0.7)
}

##########################################

png("supp_figs/fig_S11.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 1:5){
  f=families[n]
  plot_lines(Bn_imp_p[,f], Bn_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()

############################################

png("supp_figs/fig_S12.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 1:5){
  f=families[n]
  plot_lines(At_imp_p[,f], At_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()

##########################################

png("supp_figs/fig_S13.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 6:10){
  f=families[n]
  plot_lines(Bn_imp_p[,f], Bn_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()

##########################################

png("supp_figs/fig_S14.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 6:10){
  f=families[n]
  plot_lines(At_imp_p[,f], At_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()

#################################################

png("supp_figs/fig_S15.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 11:15){
  f=families[n]
  plot_lines(Bn_imp_p[,f], Bn_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()

#################################################

png("supp_figs/fig_S16.png", width=17, height=20, res=800, units="cm")
layout(matrix(c(1:14), ncol=2, byrow=T),
       heights=c(0.1,rep(1,5),0.1),
       widths=rep(2.5,2))
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Promoter")))
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels=expression(bold("Terminator")))
for (n in 11:15){
  f=families[n]
  plot_lines(At_imp_p[,f], At_imp_t[,f], dap_p[,f], dap_t[,f], families_full[n])
}
par(mar=c(0,5,0,0.5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TSS")
par(mar=c(0,0.5,0,5))
plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
text(x=0, y=0.5, labels="Position Relative to TTS")
dev.off()



######################################
######################################



#### FIGURES S18-S23 ######## Motif insertion analysis
##################################################

plot_insert <- function(delta_p_q1, delta_t_q1, delta_p_q2, delta_t_q2, delta_p_q3, delta_t_q3, family){
  x <- c(-500:499)
  y_min <- min(c(delta_p_q1, delta_t_q1))
  y_max <- max(c(delta_p_q3, delta_t_q3))
  par(mar=c(2,3,1,0.5))
  plot(x, delta_p_q2, type="l", ylim=c(y_min,y_max), axes=F, ylab="", xlab="")
  polygon(x=c(x,rev(x)), y=c(delta_p_q3,rev(delta_p_q1)), col=transparent("blue", alpha=0.3), border=NA)
  lines(x, delta_p_q2, lwd=2)
  abline(h=0, lty=2, col=transparent("grey40",alpha=0.5))
  if (max(delta_p_q2) <= 0.02){
    legend("bottomleft", pch=NA, legend=family, bty="n", text.col="grey20", cex=1.1) 
  } else {
    legend("topleft", pch=NA, legend=family, bty="n", text.col="grey20", cex=1.1)  
  }
  axis(2)
  axis(1)
  
  par(mar=c(2,0.5,1,3))
  plot(x, delta_t_q2, type="l", ylim=c(y_min,y_max), axes=F, ylab="", xlab="")
  polygon(x=c(x,rev(x)), y=c(delta_t_q3,rev(delta_t_q1)), col=transparent("blue", alpha=0.3), border=NA)
  lines(x, delta_t_q2, lwd=2)
  abline(h=0, lty=2, col=transparent("grey40",alpha=0.5))
  axis(1)
}



get_mutation_data <- function(organism, type="TFBS", quantile=0.5){
  pdf <- data.frame(position=c(-500:499))
  tdf <- data.frame(position=c(-500:499))
  for (family in families){
    df <- read.table(paste0("../result_subset/nemo/",organism,"/masked_graphpart/motif_insertion/", family, "_medium.tsv"), header=T)
    df_contr <- df[df$motif_type=="random",] 
    df <- df[df$motif_type==type,]
    df$diff <- (df$mutated_pred - df_contr$mutated_pred)
    pdiff <- c()
    tdiff <- c()
    for (p in pdf$position){
      pdiff <- c(pdiff, quantile(df$diff[df$mutated_sequence_type=="promoter" & df$relative_idx==p], probs=quantile))
      tdiff <- c(tdiff, quantile(df$diff[df$mutated_sequence_type=="terminator" & df$relative_idx==p], probs=quantile))
    }
    pdf <- cbind(pdf, pdiff)
    colnames(pdf)[ncol(pdf)] <- family
    tdf <- cbind(tdf, tdiff)
    colnames(tdf)[ncol(tdf)] <- family
  }
  for (f in families){
    R=1
    N=101 #75
    P=1
    pdf[,f] <- sgolayfilt(pdf[,f], p = P, n = N)
    tdf[,f] <- sgolayfilt(tdf[,f], p = P, n = N)
  }
  return(list(pdf,tdf))
}

make_plot <- function(filename, range){
  png(filename, width=17, height=20, res=800, units="cm")
  layout(matrix(c(1,2,3,
                  1,4,5,
                  1,6,7,
                  1,8,9,
                  1,10,11,
                  1,12,13,
                  1,14,15), ncol=3, byrow=T),
         heights=c(0.1,rep(1,5),0.1),
         widths=c(0.2,rep(2.5,2)))
  par(mar=c(0,0,0,0.0))
  plot(x=c(0,1),y=c(0,1), type="n", axes=F, ylab="", xlab="")
  text(x=0.5,y=0.5,labels="Change in Prediction w.r.t. Control (Median and Interquartile Range)", srt=90, cex=1.5)
  par(mar=c(0,3,0,0.5))
  plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
  text(x=0, y=0.5, labels=expression(bold("Promoter")))
  par(mar=c(0,0.5,0,3))
  plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
  text(x=0, y=0.5, labels=expression(bold("Terminator")))
  ###################
  for (f in families[range])
    plot_insert(pq1[,f], tq1[,f], pq2[,f], tq2[,f], pq3[,f], tq3[,f], f)
  ###################
  par(mar=c(0,5,0,0.5))
  plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
  text(x=0, y=0.5, labels="Position Relative to TSS")
  par(mar=c(0,0.5,0,5))
  plot(x=c(-500,500),y=c(0,1), type="n", axes=F, ylab="", xlab="")
  text(x=0, y=0.5, labels="Position Relative to TTS")
  dev.off()
}

organism <- "Bnapus"

q1_data <- get_mutation_data(organism, type="TFBS", quantile=0.25)
q2_data <- get_mutation_data(organism, type="TFBS", quantile=0.5)
q3_data <- get_mutation_data(organism, type="TFBS", quantile=0.75)
pq1 <- q1_data[[1]]
tq1 <- q1_data[[2]]
pq2 <- q2_data[[1]]
tq2 <- q2_data[[2]]
pq3 <- q3_data[[1]]
tq3 <- q3_data[[2]]

make_plot("supp_figs/fig_S18.png", c(1:5)) 
make_plot("supp_figs/fig_S20.png", c(6:10))  
make_plot("supp_figs/fig_S22.png", c(11:15))  


organism <- "Athaliana"

q1_data <- get_mutation_data(organism, type="TFBS", quantile=0.25)
q2_data <- get_mutation_data(organism, type="TFBS", quantile=0.5)
q3_data <- get_mutation_data(organism, type="TFBS", quantile=0.75)
pq1 <- q1_data[[1]]
tq1 <- q1_data[[2]]
pq2 <- q2_data[[1]]
tq2 <- q2_data[[2]]
pq3 <- q3_data[[1]]
tq3 <- q3_data[[2]]

make_plot("supp_figs/fig_S19.png", c(1:5)) 
make_plot("supp_figs/fig_S21.png", c(6:10))  
make_plot("supp_figs/fig_S23.png", c(11:15))  
  
