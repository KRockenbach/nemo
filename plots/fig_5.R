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


# processing steps:
########################
# importance:
# 1) importance for each gene is calculated as the sum of absolute attributions across bases at each position
# 2) global mean and sd across all genes is calculated at each position (taking masked segments into account)
# 3) mean across each group of genes is calculated (taking masked segments into account)
# 4) savitzky-golay filter is applied to global mean, global sd, and group means
# 5) global mean is mean-normalized (mean across both sequences is subtracted from each position)
# 6) group means are mean-normalized (mean across both sequences is subtracted from each position)
# 7) normalized global mean is subtracted from normalized group means
# 8) difference is divided by global standard deviation
# 9) color values are scaled globally for eahc sequence 
    #(min-max scaling to the range [0,1], 
    #where the 0 represents the minumum across groups and positions 
    #and 1 represents the maximum across groups and positions)


# DAP-seq:
# Min-max scaling to range [0,1] within each group (by row) across both sequences! 



organisms <- c("Bnapus", "Athaliana")

families <- c("MYBrelated", "Homeobox", "C2C2dof", "WRKY",
              "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
              "NAC", "G2like", "HSF", "AP2EREBP", "C2C2gata")
families_full <- c("MYB-related", "Homeobox", "C2C2-Dof", "WRKY",
                   "TCP", "bZIP", "MYB", "MADS", "bHLH", "Trihelix",
                   "NAC", "G2-like", "HSF", "AP2/EREBP", "C2C2-GATA")


fig_lab_cex=1.5


plot_heatmap <- function(pdf, tdf, filter=FALSE, palette="custom",
                         descriptor="", norm_type="positive", legend_labs=c("Low", "High"), 
                         label=expression(bold("C")), byrow=T, organism){ #promoter DF, terminator Df
  par(mar=c(2.5,0,1.5,0), tck=-0.03, mgp=c(1.5,0.4,0))
  plot(x=c(-460:460), y=c(0:-920), type="n", 
       axes=F, xlab="", ylab="", ylim=c(-1500,0))
  fig_label(label, cex=fig_lab_cex)
  # joint normalization and smoothing of promoter and terminator data
  for (f in families){
    if(filter){
      pdf[,f] <- sgolayfilt(pdf[,f], p = 1, n = 45) #75)
      tdf[,f] <- sgolayfilt(tdf[,f], p = 1, n = 45) #75) 
    }
    if(norm_type=="center"){
      norm.cat <- min_max.zero_center.scale(c(pdf[,f], tdf[,f]))
      pdf[,f] <- norm.cat[1:1001]
      tdf[,f] <- norm.cat[1002:2002] 
    } else if(norm_type=="positive"){
      norm.cat <- min_max.positive.scale(c(pdf[,f], tdf[,f]))
      pdf[,f] <- norm.cat[1:1001]
      tdf[,f] <- norm.cat[1002:2002]
    }
  }  
  if (norm_type=="positive_global"){
    norm.cat <- min_max.positive.scale(as.matrix(rbind(pdf[,families], tdf[,families])))
    pdf[,families] <- norm.cat[1:1001,]
    tdf[,families] <- norm.cat[1002:2002,]
  } else if(norm_type=="center_global"){
    norm.cat <- min_max.zero_center.scale(as.matrix(rbind(pdf[,families], tdf[,families])))
    pdf[,families] <- norm.cat[1:1001,]
    tdf[,families] <- norm.cat[1002:2002,] 
  }
  for (f in 1:length(families)){
    y1 <- (-1)*(f-1)*100
    y2 <- (-1)*(f)*100
    text(x=300, y=(y1+y2)/2, labels=families_full[f], adj=1, cex=0.75)
  }
  
  if (organism=="Athaliana"){
    text(x=300, y=120, labels=expression(italic("A. thaliana")), adj=1, cex=1, col="grey40")
  } else {
    text(x=300, y=120, labels=expression(italic("B. napus")), adj=1, cex=1, col="grey40")
  }
  
  for (seq in c("Promoter", "Terminator")){
    if (seq == "Promoter"){
      df <- pdf
      ref_point <- "TSS"
    } else {
      df <- tdf
      ref_point <- "TTS"
    }
    if (seq == "Promoter"){
      par(mar=c(2.5,0.1,1.5,0.3), cex.axis=0.6, xpd=T)
    } else {
      par(mar=c(2.5,0.3,1.5,0.1), cex.axis=0.6, xpd=T)
    }
    plot(x=c(-460:460), y=c(0:-920), type="n", 
         axes=F, xlab="", ylab="", ylim=c(-1500,0), cex.lab=lab_cex)
    if (byrow){ # scale color values by row (column in df; between positions)
      for (f in 1:15){
        fam <- families[f]
        y1 <- (-1)*(f-1)*100
        y2 <- (-1)*f*100
        cols <- get_colors(df[,fam], palette=palette)
        for (idx in 1:1001){
          x1 <- (-500) + (idx-1)
          x2 <- (-500) + idx
          polygon(x=c(x1,x1,x2,x2), y=c(y1,y2,y2,y1), border=NA, col=cols[idx]) 
        }
        polygon(x=c(-500,-500,500,500), y=c(y1,y2,y2,y1), lwd=1)
      } 
    } else { # scale color values globally
      min <- min(df[,families])
      max <- max(df[,families])
      for (idx in 1:1001){
        x1 <- (-500) + (idx-1)
        x2 <- (-500) + idx
        cols <- get_colors(as.matrix(df[idx,families]), palette=palette)
        for (f in 1:15){
          y1 <- (-1)*(f-1)*100
          y2 <- (-1)*f*100
          polygon(x=c(x1,x1,x2,x2), y=c(y1,y2,y2,y1), border=NA, col=cols[f]) 
        }
      }
      for (f in 1:15){
        y1 <- (-1)*(f-1)*100
        y2 <- (-1)*f*100
        polygon(x=c(-500,-500,500,500), y=c(y1,y2,y2,y1), lwd=1)
      }
    }
    lines(x=c(0,0), y=c(-1500,0), col=transparent("grey20",0.6), lwd=2)
    axis(side=1, at=seq(-500,500,50), las=2, cex.axis=ax_cex)
    title(xlab=paste0("Position Relative to ", ref_point), cex.lab=lab_cex)
    title(main=seq, line=0.02, cex.main=main_cex)
  }
  par(mar=c(1.5,0.5,1.5,0.5))
  plot(x=c(-60:60), y=c(0:120), type="n", 
       axes=F, xlab="", ylab="", ylim=c(0,120), cex.lab=lab_cex)
  if (palette=="custom"){
    cols <- get_colors(seq(-1,1,2/999), palette="custom") 
  } else {
    cols <- get_colors(seq(0,1,1/999), palette=palette)
  }
  #legend width
  lw <- 15
  #legend height
  lh <- 0.05
  for (i in 1:1000){
    polygon(x=c(-1,-1,1,1)*lw, y=c(25+((i-1)*lh), 25+(i*lh), 25+(i*lh), 25+((i-1)*lh)), col=cols[i], border=NA)
  }
  lines(x=c(-lw-1,lw+1), y=c(25,25))
  lines(x=c(-lw-1,lw+1), y=c(75,75))
  text(x=c(0,0), y=c(20,80), labels=legend_labs, cex=0.6, adj=0.5, col="grey40")
  par(xpd=T)
  text(x=0, y=110, labels=descriptor, cex=0.6)
}



png("fig_4.png", height=20, width=17, res=1200, units="cm")
par(xpd=T)
layout(matrix(c(1,2,3,4,
                5,6,7,8,
                9,10,11,12,
                13,14,15,16,
                17,18,19,20), ncol=4, byrow=T),
       heights=rep(1,5),
       widths=c(1.25,3.5,3.5,0.75))


seqs <- c("promoter", "terminator")


lab_cex=0.75
ax_cex=0.6
main_cex=0.85

pdf <- data.table::transpose(read.table("../data/plot_data/DAPseq/DAPSeq.TF.promoter.tsv", header=F), make.names = 1)
pdf <- pdf[,families]
position <- as.numeric(rownames(pdf))
pdf <- cbind(position, pdf)
tdf <- data.table::transpose(read.table("../data/plot_data/DAPseq/DAPSeq.TF.terminator.tsv", header=F), make.names = 1)
tdf <- tdf[,families]
tdf <- cbind(position, tdf)

plot_heatmap(pdf,tdf, palette="inferno", descriptor="DAP-Seq\nSignal\n(Scaled by Row)", 
             label=expression(bold("A")), norm_type="positive", organism="Athaliana")

for (o in 1:2){
    organism <- organisms[o]
    
    
    pdf <- read.table(paste0("../results/nemo/",organism,
                             "/masked_graphpart/attribs/promoter_TF_importance.tsv"),
                      header=T)
    tdf <- read.table(paste0("../results/nemo/",organism,
                             "/masked_graphpart/attribs/terminator_TF_importance.tsv"),
                      header=T)
    tdf$position <- tdf$position + 1 # shift position, because TTS is part of upstream sequence
    # trim to +- 500 bp
    pdf <- pdf[501:1501,]
    tdf <- tdf[501:1501,]
    pdf <- pdf[,c("position","mean", "sd", families)]
    tdf <- tdf[,c("position","mean", "sd", families)]
    
    
    # savitzky-golay filtering
    for (c in 2:ncol(pdf)){
      R=1
      N=45 #75
      P=1
      pdf[,c] <- sgolayfilt(pdf[,c], p = P, n = N)
      tdf[,c] <- sgolayfilt(tdf[,c], p = P, n = N)
    }
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

    if (organism=="Athaliana"){
      label=expression(bold("C"))
    } else {
      label=expression(bold("B"))
    }

    plot_heatmap(pdf,tdf, filter=FALSE, palette="inferno", norm_type="positive_global",
                 descriptor="Normalized\nImportance\n(Scaled Globally)",
                 byrow=F, label=label, organism=organism)
}



for (organism in c("Bnapus","Athaliana")){   
    pdf <- read.table(paste0("../data/plot_data/", organism, "_promoter_TF_inserted.tsv"), header=T)
    tdf <- read.table(paste0("../data/plot_data/", organism, "_terminator_TF_inserted.tsv"), header=T)
    tdf$x <- tdf$x + 1 # shift position, because TTS is part of upstream sequence
    # trim to +- 500 bp
    pdf <- pdf[pdf$x %in% c(-500:500),]
    tdf <- tdf[tdf$x %in% c(-500:500),]
    
    if (organism=="Athaliana"){
      label=expression(bold("E"))
    } else {
      label=expression(bold("D"))
    }
    
    plot_heatmap(pdf,tdf, filter=TRUE, palette="custom", 
                 descriptor="Change in\nExpression\n(Scaled Globally)", 
                 legend_labs=c("Negative","Positive"), label=label, norm_type="center_global", organism=organism)

}

dev.off()
